# -----------------------------
# QC Filter UI
# -----------------------------
qcFilterUI <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Quality Control Filtering"),
    
    # Intro
    p("This module allows you to filter low-quality cells for each sample based on standard QC metrics and doublet predictions."),
    hr(),
    
    # QC metrics explanation
    h3("QC Metrics Explained"),
    tags$ul(
      tags$li(strong("nFeature_RNA:"), " Number of detected genes per cell. Low = empty droplets, high = multiplets."),
      tags$li(strong("nCount_RNA:"), " Total RNA counts per cell."),
      tags$li(strong("percent.mt:"), " Percentage of mitochondrial reads. High values may indicate stressed/dying cells."),
      tags$li(strong("percent.ribo:"), " Percentage of ribosomal reads. High values may bias clustering."),
      tags$li(strong("doublet_predictions:"), " Predicted doublets. Removing them improves downstream analysis.")
    ),
    hr(),
    
    # Filtering guidance
    h3("How to Set Filters"),
    tags$ul(
      tags$li("Use violin plots to visualize distribution of each metric."),
      tags$li("Exclude extreme outliers."),
      tags$li("Optionally remove predicted doublets.")
    ),
    p("After adjusting thresholds, click ", strong("Apply Filter"), " for each sample."),
    hr(),
    
    # Step 1: per-sample inspection
    h3("Step 1: Inspect QC per Sample"),
    p("Adjust thresholds and doublet removal per sample."),
    uiOutput(ns("per_sample_panels")),
    hr(),
    
    # Step 2: Summary
    h3("Step 2: Review Summary"),
    DT::DTOutput(ns("qc_summary_table")),
    br(),
    verbatimTextOutput(ns("filter_info"))
  )
}

# -----------------------------
# QC Filter Server
# -----------------------------
qcFilterServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
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
    
    # ------------------ Ensure QC and doublet info ------------------
    observeEvent(app_data$seurat_obj, {
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      
      # percent.mt
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        patt <- .guess_mito_pattern(obj)
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = patt)
      }
      
      # percent.ribo
      if (!"percent.ribo" %in% colnames(obj@meta.data)) {
        ribo_genes <- grep("^RPS|^RPL|^Rps|^Rpl", rownames(obj), value = TRUE)
        if (length(ribo_genes) > 0) {
          obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_genes)
        } else {
          warning("⚠️ No ribosomal genes found — setting percent.ribo = 0") # ⚠️ FIX
          obj$percent.ribo <- 0
        }
      }
      
      app_data$seurat_obj <- obj
    }, ignoreInit = FALSE)
    
    # ------------------ Per-sample UI panels ------------------
    output$per_sample_panels <- renderUI({
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      grp <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = grp)
      if (!length(parts)) return(tags$em("No samples found."))
      
      tagList(lapply(names(parts), function(s) {
        sid <- .safe_id(s)
        so <- parts[[s]]
        
        # Add missing QC metrics if necessary
        if (!"percent.mt" %in% colnames(so@meta.data)) so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = .guess_mito_pattern(so))
        if (!"percent.ribo" %in% colnames(so@meta.data)) so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^RPS|^RPL")
        
        # Default thresholds
        qf_lo <- suppressWarnings(quantile(so$nFeature_RNA, 0.05, na.rm = TRUE)) %||% 200
        qf_hi <- suppressWarnings(quantile(so$nFeature_RNA, 0.95, na.rm = TRUE)) %||% 2500
        qmt   <- suppressWarnings(quantile(so$percent.mt, 0.95, na.rm = TRUE)) %||% 5
        qribo <- suppressWarnings(quantile(so$percent.ribo, 0.95, na.rm = TRUE)) %||% 50
        
        local({
          .sid <- sid; .so <- so; .s <- s
          output[[paste0("vln_", .sid)]] <- renderPlot({
            p1 <- VlnPlot(.so, features = "nFeature_RNA", pt.size = 0.2) + ggtitle("nFeature_RNA")
            p2 <- VlnPlot(.so, features = "nCount_RNA", pt.size = 0.2) + ggtitle("nCount_RNA")
            p3 <- VlnPlot(.so, features = "percent.mt", pt.size = 0.2) + ggtitle("percent.mt")
            p4 <- VlnPlot(.so, features = "percent.ribo", pt.size = 0.2) + ggtitle("percent.ribo")
            (p1 | p2 | p3 | p4) + plot_layout(ncol = 4)
          }, height = 600)
        })
        
        wellPanel(
          h4(sprintf("Sample: %s", s)),
          plotOutput(session$ns(paste0("vln_", sid)), height = "600px"),
          fluidRow(
            column(3, numericInput(session$ns(paste0("minF_", sid)), "Min #features", round(qf_lo), min = 0)),
            column(3, numericInput(session$ns(paste0("maxF_", sid)), "Max #features", round(qf_hi), min = 0)),
            column(3, numericInput(session$ns(paste0("maxMT_", sid)), "Max % mito", round(qmt,1), min = 0, max = 100, step = 0.1)),
            column(3, numericInput(session$ns(paste0("maxRibo_", sid)), "Max % ribo", round(qribo,1), min = 0, max = 100, step = 0.1))
          ),
          checkboxInput(session$ns(paste0("removeDoublets_", sid)), "Remove predicted doublets", value = TRUE),
          actionButton(session$ns(paste0("apply_", sid)), "Apply Filter for This Sample", class = "btn-primary")
        )
      }))
    })
    
    # ------------------ Apply filters ------------------
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
            
            if (is.null(s_name) || !(s_name %in% names(parts_now))) {
              showNotification("❌ Could not resolve sample.", type="error")
              return()
            }
            
            # thresholds
            minF <- input[[paste0("minF_", .sid)]] %||% 200
            maxF <- input[[paste0("maxF_", .sid)]] %||% 2500
            maxMT <- input[[paste0("maxMT_", .sid)]] %||% 5
            maxRibo <- input[[paste0("maxRibo_", .sid)]] %||% 50
            removeDoublets <- input[[paste0("removeDoublets_", .sid)]] %||% TRUE
            
            so <- parts_now[[s_name]]
            
            # Run doublet detection only if removeDoublets is TRUE and column doesn't exist
            if (removeDoublets && !"doublet_predictions" %in% colnames(so@meta.data)) {
              sce <- as.SingleCellExperiment(so)
              sce <- scDblFinder(sce)
              so$doublet_scores <- sce$scDblFinder.score
              so$doublet_predictions <- sce$scDblFinder.class
            }
            # ⚠️ FIX: Only apply doublet filter if column exists; else keep all TRUE
            if ("doublet_predictions" %in% colnames(so@meta.data) && removeDoublets) {
              doublet_condition <- so$doublet_predictions == "singlet"
              doublet_condition[is.na(doublet_condition)] <- FALSE   # ⚠️ FIX: replace NA with FALSE
            } else {
              doublet_condition <- rep(TRUE, ncol(so))              # ⚠️ FIX: no doublet column → all pass
            }
            
            # ⚠️ FIX: ensure other metrics don’t have NA
            keep_features <- !is.na(so$nFeature_RNA) & so$nFeature_RNA > minF & so$nFeature_RNA < maxF
            keep_mt       <- !is.na(so$percent.mt) & so$percent.mt < maxMT
            keep_ribo     <- !is.na(so$percent.ribo) & so$percent.ribo < maxRibo
            keep <- keep_features & keep_mt & keep_ribo & doublet_condition
            
            filtered <- subset(so, cells = colnames(so)[keep])  # ⚠️ FIX: subset with logical index
            
            parts_now[[s_name]] <- filtered
            nonempty <- sapply(parts_now, function(x) ncol(x)) > 0
            parts_now <- parts_now[nonempty]
            
            # --- NEW: check if all cells were filtered out ---
            if (length(parts_now) == 0) {
              showNotification("⚠️ All cells were filtered out. Please relax your thresholds.", type = "warning")
              return()
            }
            
            kept_names <- names(parts_now)
            
            merged <- if (length(parts_now) == 1) {
              one <- parts_now[[1]]
              if (ncol(one) > 0) {
                colnames(one) <- paste0(kept_names[1], "_", colnames(one))
              }
              one
            } else {
              merge(
                x = parts_now[[1]],
                y = parts_now[-1],
                add.cell.ids = kept_names,
                project = "QC_filtered"
              )
            }
            
            app_data$seurat_obj <- merged
            
            
            # Update summary table
            cur_parts <- SplitObject(merged, split.by = "orig.ident")
            summary_df <- data.frame(
              sample = names(cur_parts),
              cells = sapply(cur_parts, ncol),
              stringsAsFactors = FALSE
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
    
    output$qc_summary_table <- DT::renderDT({ NULL })
    output$filter_info <- renderText({ "" })
    
  })
}
