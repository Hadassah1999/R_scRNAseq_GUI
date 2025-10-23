annotationUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Cell Annotation",
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          ns("reference_choice"), "Reference Genome:",
          choices = c("Mouse" = "mouse", "Human" = "human"),
          selected = "mouse"
        ),
        actionButton(ns("annotate_btn"), "Run SingleR Annotation"),
        uiOutput(ns("dataset_filter_ui")),  # <- added filter UI
        checkboxInput(ns("show_labels"), "Show annotation labels on UMAP", value = TRUE),
        downloadButton(ns("download_plot"), "Download DimPlot")
      ),
      mainPanel(
        plotOutput(ns("plot")),
        plotOutput(ns("pie_plot")),
        tableOutput(ns("meta_preview"))
      )
    )
  )
}

annotationServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    running <- reactiveVal(FALSE)  # <-- re-entry guard
    
    # --- helpers ---
    .pick_assay <- function(obj) {
      assays_slot <- tryCatch(names(obj@assays), error = function(e) character(0))
      assays_func <- tryCatch(names(Seurat::Assays(obj)), error = function(e) character(0))
      assays_all  <- unique(c(assays_slot, assays_func))
      if (length(assays_all) == 0) stop("This Seurat object has no assays.")
      da <- suppressWarnings(tryCatch(Seurat::DefaultAssay(obj), error = function(e) NA_character_))
      if (!is.null(da) && nzchar(da) && da %in% assays_all) return(da)
      Seurat::DefaultAssay(obj) <- assays_all[1]
      assays_all[1]
    }
    
    .layers <- function(obj, assay) {
      tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
    }
    .has_data_layer <- function(obj, assay) {
      any(grepl("^data(\\.|$)", .layers(obj, assay)))
    }
    
    # Collect ALL "data" / "data.*" layers and sparse-cbind them (no JoinLayers)
    .extract_log_mat <- function(obj, assay_for_join) {
      if (!.has_data_layer(obj, assay_for_join)) {
        obj <- Seurat::NormalizeData(obj, assay = assay_for_join, verbose = FALSE)
      }
      lyr_names <- .layers(obj, assay_for_join)
      data_like <- lyr_names[grepl("^data(\\.|$)", lyr_names)]
      if (length(data_like) == 0) stop("No 'data' layers found after normalization.")
      
      mats <- lapply(data_like, function(ln) {
        m <- SeuratObject::LayerData(obj[[assay_for_join]], layer = ln)
        if (!inherits(m, "dgCMatrix")) m <- methods::as(m, "dgCMatrix")
        m
      })
      mat <- do.call(cbind, mats)  # still sparse
      
      # align to object’s cell order
      cells_all <- colnames(obj)
      if (is.null(rownames(mat))) stop("'data' matrix has no gene rownames.")
      if (is.null(colnames(mat))) stop("'data' matrix has no cell colnames.")
      keep <- intersect(cells_all, colnames(mat))
      mat  <- mat[, keep, drop = FALSE]
      miss <- setdiff(cells_all, colnames(mat))
      if (length(miss)) {
        add <- Matrix::Matrix(0, nrow = nrow(mat), ncol = length(miss), sparse = TRUE)
        rownames(add) <- rownames(mat); colnames(add) <- miss
        mat <- cbind(mat, add)
      }
      mat[, cells_all, drop = FALSE]
    }
    
    .choose_group_field <- function(obj) {
      if ("orig.ident" %in% colnames(obj@meta.data)) return("orig.ident")
      if ("sample"     %in% colnames(obj@meta.data)) return("sample")
      stop("No 'orig.ident' or 'sample' in metadata for per-sample fallback.")
    }
    
    # Keep test sparse; only ref (small) is dense
    .run_singler_matrix <- function(test_log, ref_se) {
      ref_log <- SummarizedExperiment::assay(ref_se, "logcounts")
      common  <- intersect(rownames(test_log), rownames(ref_log))
      if (length(common) < 50) stop(sprintf("Too few overlapping genes with reference: %d", length(common)))
      
      test_c <- test_log[common, , drop = FALSE]                 # dgCMatrix (sparse)
      ref_c  <- as.matrix(ref_log[common, , drop = FALSE])       # dense OK (small)
      lab_c  <- ref_se$label.main
      
      SingleR::SingleR(
        test  = test_c,
        ref   = ref_c,
        labels= lab_c,
        assay.type.test = NULL,
        assay.type.ref  = NULL
      )
    }
    
    # --- SingleR annotation ---
    observeEvent(input$annotate_btn, {
      if (isTRUE(running())) return(invisible(NULL))  # guard against re-entry
      running(TRUE)
      on.exit(running(FALSE), add = TRUE)
      
      obj <- app_data$dataset %||% app_data$seurat_obj
      req(obj)
      showNotification("🔄 Running SingleR…", type = "message", id = "anno_progress", duration = NULL)
      
      tryCatch({
        aa <- .pick_assay(obj)
        if (!.has_data_layer(obj, aa)) obj <- Seurat::NormalizeData(obj, assay = aa, verbose = FALSE)
        log_all <- .extract_log_mat(obj, assay_for_join = aa)   # sparse & assembled from data.*
        
        ref_se <- switch(input$reference_choice,
                         "mouse" = celldex::MouseRNAseqData(),
                         "human" = celldex::HumanPrimaryCellAtlasData())
        
        pred <- tryCatch(.run_singler_matrix(log_all, ref_se), error = function(e) NULL)
        if (!is.null(pred)) {
          obj$cell_ann <- pred$labels
        } else {
          # fallback per-group (still sparse)
          grp <- .choose_group_field(obj)
          labels <- rep(NA_character_, ncol(obj)); names(labels) <- colnames(obj)
          for (g in unique(obj@meta.data[[grp]])) {
            idx <- which(obj@meta.data[[grp]] == g)
            if (!length(idx)) next
            mat_g <- log_all[, idx, drop = FALSE]
            pred_g <- tryCatch(.run_singler_matrix(mat_g, ref_se), error = function(e) NULL)
            if (!is.null(pred_g)) labels[idx] <- pred_g$labels
          }
          obj$cell_ann <- labels
        }
        
        app_data$dataset <- obj
        app_data$seurat_obj <- obj
        
        
        removeNotification("anno_progress")
        showNotification("✅ Annotation done.", type = "message")
      }, error = function(e) {
        removeNotification("anno_progress")
        showNotification(paste("❌ Annotation error:", e$message), type = "error")
      })
    }, ignoreInit = TRUE)  # <-- don’t fire on app load
    
    # --- Dataset selection filter ---
    output$dataset_filter_ui <- renderUI({
      obj <- app_data$dataset %||% app_data$seurat_obj
      req(obj); req("cell_ann" %in% colnames(obj@meta.data))
      grp_col <- if ("orig.ident" %in% colnames(obj@meta.data)) "orig.ident" else "sample"
      choices <- unique(obj@meta.data[[grp_col]])
      
      tagList(
        radioButtons(session$ns("plot_type"), "Choose Plot Type:",
                     choices = c("All datasets combined" = "all", "Select datasets" = "subset"),
                     selected = "all"),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'subset'", session$ns("plot_type")),
          checkboxGroupInput(session$ns("dataset_choice"), "Select dataset(s):",
                             choices = choices, selected = choices[1])
        )
      )
    })
    
    # --- Reactive filtered object ---
    filtered_obj <- reactive({
      obj <- app_data$dataset %||% app_data$seurat_obj
      req(obj); req("cell_ann" %in% colnames(obj@meta.data))
      grp_col <- if ("orig.ident" %in% colnames(obj@meta.data)) "orig.ident" else "sample"
      
      if (!is.null(input$plot_type) && input$plot_type == "subset" && !is.null(input$dataset_choice)) {
        cells_to_keep <- rownames(obj@meta.data[obj@meta.data[[grp_col]] %in% input$dataset_choice, ])
        obj <- subset(obj, cells = cells_to_keep)
      }
      obj
    })
    
    # --- Plots ---
    output$plot <- renderPlot({
      obj <- filtered_obj()
      req(obj)
      DimPlot(obj, group.by = "cell_ann", label = input$show_labels, repel = TRUE) +
        ggtitle("SingleR Annotation")
    })
    
    output$pie_plot <- renderPlot({
      obj <- filtered_obj()
      req(obj)
      df <- as.data.frame(table(obj$cell_ann))
      colnames(df) <- c("CellType", "Count")
      df$Percent <- round(100 * df$Count / sum(df$Count), 1)
      df$LegendLabel <- paste0(df$CellType, " (", df$Percent, "%)")
      df$LegendLabel <- factor(df$LegendLabel, levels = df$LegendLabel[order(-df$Percent)])
      
      ggplot(df, aes(x = "", y = Count, fill = LegendLabel)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() +
        ggtitle("Cell Type Distribution (SingleR)") +
        theme(legend.title = element_blank())
    })
    
    # --- Download DimPlot ---
    output$download_plot <- downloadHandler(
      filename = function() paste0("cell_annotation_", Sys.Date(), ".png"),
      content = function(file) {
        obj <- filtered_obj()
        req(obj)
        png(file, width = 1200, height = 900)
        print(DimPlot(obj, group.by = "cell_ann", label = TRUE, repel = TRUE))
        dev.off()
      }
    )
  })
}