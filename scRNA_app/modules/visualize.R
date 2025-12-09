qcFilterUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    h2("Quality Control Filtering"),
    
    p("This module allows you to filter low-quality cells for each sample based on standard QC metrics and doublet predictions."),
    hr(),
    
    h3("QC Metrics Explained"),
    tags$ul(
      tags$li(strong("nFeature_RNA:"), " Number of detected genes per cell. Low = empty droplets, high = multiplets."),
      tags$li(strong("nCount_RNA:"), " Total RNA counts per cell."),
      tags$li(strong("percent.mt:"), " Percentage of mitochondrial reads. High values may indicate stressed/dying cells."),
      tags$li(strong("percent.ribo:"), " Percentage of ribosomal reads. High values may bias clustering."),
      tags$li(strong("doublet_predictions:"), " Predicted doublets (computed automatically when data is loaded).")
    ),
    hr(),
    
    # Remove doublets button
    actionButton(ns("remove_doublets_global"), "Remove Doublets", class = "btn-warning"),
    hr(),
    
    h3("How to Set Filters"),
    tags$ul(
      tags$li("Use violin plots to visualize distribution of each metric."),
      tags$li("Exclude extreme outliers.")
    ),
    p("After adjusting thresholds, click ", strong("Apply Filter"), " for each sample."),
    hr(),
    
    div(id = ns("qc_panel"),
          h3("Step 1: Inspect QC per Sample"),
          p("Adjust thresholds per sample."),
          uiOutput(ns("per_sample_panels")),
          hr(),
          h3("Step 2: Review Summary"),
          DT::DTOutput(ns("qc_summary_table")),
          br(),
          verbatimTextOutput(ns("filter_info"))
      )
    )
}

qcFilterServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
    .safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
    .choose_group_col <- function(obj) {
      if ("orig.ident" %in% colnames(obj@meta.data)) return("orig.ident")
      if ("sample" %in% colnames(obj@meta.data)) return("sample")
      stop("No 'orig.ident' or 'sample' column found in metadata.")
    }
    .guess_mito_pattern <- function(obj) {
      feats <- rownames(obj)
      if (length(feats) && any(grepl("^MT-", feats))) "^MT-" else "^mt-"
    }
    
    # Ensure QC metrics are calculated immediately
    observe({
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        patt <- .guess_mito_pattern(obj)
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = patt)
      }
      if (!"percent.ribo" %in% colnames(obj@meta.data)) {
        ribo_genes <- grep("^RPS|^RPL|^Rps|^Rpl", rownames(obj), value = TRUE)
        obj[["percent.ribo"]] <- if (length(ribo_genes)) PercentageFeatureSet(obj, features = ribo_genes) else 0
      }
      app_data$seurat_obj <- obj
    })
    
    # Render per-sample panels (violin plots + thresholds)
    output$per_sample_panels <- renderUI({
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      grp <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = grp)
      if (!length(parts)) return(tags$em("No samples found."))
      
      tagList(lapply(names(parts), function(s) {
        sid <- .safe_id(s)
        so <- parts[[s]]
        
        # defaults
        qf_lo <- suppressWarnings(quantile(so$nFeature_RNA, 0.05, na.rm=TRUE)) %||% 200
        qf_hi <- suppressWarnings(quantile(so$nFeature_RNA, 0.95, na.rm=TRUE)) %||% 2500
        qmt   <- suppressWarnings(quantile(so$percent.mt, 0.95, na.rm=TRUE)) %||% 5
        qribo <- suppressWarnings(quantile(so$percent.ribo, 0.95, na.rm=TRUE)) %||% 50
        
        local({
          .sid <- sid; .so <- so; .s <- s
          output[[paste0("vln_", .sid)]] <- renderPlot({
            p1 <- VlnPlot(.so, features="nFeature_RNA", pt.size=0.2)+ggtitle("nFeature_RNA")
            p2 <- VlnPlot(.so, features="nCount_RNA", pt.size=0.2)+ggtitle("nCount_RNA")
            p3 <- VlnPlot(.so, features="percent.mt", pt.size=0.2)+ggtitle("percent.mt")
            p4 <- VlnPlot(.so, features="percent.ribo", pt.size=0.2)+ggtitle("percent.ribo")
            (p1|p2|p3|p4)+plot_layout(ncol=4)
          }, height=600)
        })
        
        wellPanel(
          h4(sprintf("Sample: %s", s)),
          plotOutput(session$ns(paste0("vln_", sid)), height="600px"),
          fluidRow(
            column(3, numericInput(session$ns(paste0("minF_", sid)), "Min #features", round(qf_lo), min=0)),
            column(3, numericInput(session$ns(paste0("maxF_", sid)), "Max #features", round(qf_hi), min=0)),
            column(3, numericInput(session$ns(paste0("maxMT_", sid)), "Max % mito", round(qmt,1), min=0, max=100, step=0.1)),
            column(3, numericInput(session$ns(paste0("maxRibo_", sid)), "Max % ribo", round(qribo,1), min=0, max=100, step=0.1))
          ),
          actionButton(session$ns(paste0("apply_", sid)), "Apply Filter for This Sample", class="btn-primary")
        )
      }))
    })
    
    # Apply per-sample filters (thresholds)
    registered <- reactiveVal(character())
    observe({
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      grp <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = grp)
      if (!length(parts)) return()
      
      sids <- sapply(names(parts), .safe_id)
      todo <- setdiff(sids, registered())
      if (!length(todo)) return()
      sid2name <- setNames(names(parts), sids)
      
      for (sid in todo) {
        btn_id <- paste0("apply_", sid)
        local({
          .sid <- sid
          observeEvent(input[[btn_id]], {
            req(app_data$seurat_obj)
            obj_now <- app_data$seurat_obj
            parts_now <- SplitObject(obj_now, split.by = .choose_group_col(obj_now))
            s_name <- sid2name[[.sid]]
            if (is.null(s_name) || !(s_name %in% names(parts_now))) return()
            
            so <- parts_now[[s_name]]
            
            minF <- input[[paste0("minF_", .sid)]] %||% 200
            maxF <- input[[paste0("maxF_", .sid)]] %||% 2500
            maxMT <- input[[paste0("maxMT_", .sid)]] %||% 5
            maxRibo <- input[[paste0("maxRibo_", .sid)]] %||% 50
            
            keep <- !is.na(so$nFeature_RNA) & so$nFeature_RNA>minF & so$nFeature_RNA<maxF &
              !is.na(so$percent.mt) & so$percent.mt<maxMT &
              !is.na(so$percent.ribo) & so$percent.ribo<maxRibo
            
            n_before <- ncol(so)
            n_after  <- sum(keep)
            n_removed <- n_before - n_after
            
            showNotification(
              paste0("✓ Filter applied to sample '", s_name, 
                     "'. Removed ", n_removed, 
                     " low-quality cells (", n_after, " kept)."),
              type = "message",
              duration = 6
            )
            
            parts_now[[s_name]] <- subset(so, cells=colnames(so)[keep])
            
            # merge
            shared_cols2 <- Reduce(intersect, lapply(parts_now, function(x) colnames(x@meta.data)))
            parts_now <- lapply(parts_now, function(x) { x@meta.data <- x@meta.data[, shared_cols2, drop=FALSE]; x })
            merged <- Reduce(function(a,b) merge(a,b), parts_now)
            app_data$seurat_obj <- merged
            
            # update summary table
            cur_parts <- SplitObject(merged, split.by = .choose_group_col(merged))
            summary_df <- data.frame(
              sample = names(cur_parts),
              cells = sapply(cur_parts, ncol),
              stringsAsFactors=FALSE
            )
            output$qc_summary_table <- DT::renderDT({ 
              DT::datatable(summary_df, options=list(pageLength=10), rownames=FALSE) 
            })
            output$filter_info <- renderText({ paste("Total cells after filtering:", ncol(merged)) })
            
          }, ignoreInit=TRUE)
        })
      }
      registered(c(registered(), todo))
    })
    
    # initial empty outputs
    output$qc_summary_table <- DT::renderDT({ NULL })
    output$filter_info <- renderText({ "" })
    
    # Remove doublets button still works
    observeEvent(input$remove_doublets_global, {
      req(app_data$seurat_obj)
      notif_id <- showNotification("Computing and removing doublets for all samples…",
                                   duration = NULL, type = "message")
      
      obj <- tryCatch(RenameCells(app_data$seurat_obj, add.cell.id = TRUE),
                      error = function(e) app_data$seurat_obj)
      split_col <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = split_col)
      
      removed_counts <- list()
      for (nm in names(parts)) {
        so <- parts[[nm]]
        sce <- as.SingleCellExperiment(so)
        sce <- scDblFinder(sce)
        so$doublet_scores <- sce$scDblFinder.score
        so$doublet_predictions <- sce$scDblFinder.class
        
        n_doublets <- sum(so$doublet_predictions == "doublet", na.rm = TRUE)
        removed_counts[[nm]] <- n_doublets
        
        keep <- so$doublet_predictions == "singlet"
        keep[is.na(keep)] <- FALSE
        parts[[nm]] <- subset(so, cells = colnames(so)[keep])
      }
      
      # Harmonize and merge
      shared_cols <- Reduce(intersect, lapply(parts, function(x) colnames(x@meta.data)))
      parts <- lapply(parts, function(x) { x@meta.data <- x@meta.data[, shared_cols, drop=FALSE]; x })
      app_data$seurat_obj <- Reduce(function(a,b) merge(a,b), parts)
      
      removeNotification(notif_id)
      total_removed <- sum(unlist(removed_counts))
      message_text <- paste0(
        "✓ Doublets removed: ", total_removed,
        "\n", paste0(names(removed_counts), ": ", unlist(removed_counts), collapse = " | ")
      )
      
      showNotification(message_text, type = "message", duration = 8)
    })
    
  })
}
