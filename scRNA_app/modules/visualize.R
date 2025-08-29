qcFilterUI <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Quality Control Filtering"),
    
    # Intro section
    p("This module allows you to filter low-quality cells for each sample. The filtering is based on standard QC metrics."),
    hr(),
    
    # Explanation of metrics
    h3("QC Metrics Explained"),
    tags$ul(
      tags$li(
        strong("nFeature_RNA:"), 
        " Number of detected genes per cell. Keep the main cluster of cells and exclude extreme low/high values. 
         Low values may represent empty droplets, very high values may represent multiplets."
      ),
      tags$li(
        strong("nCount_RNA:"), 
        " Total RNA molecule counts per cell. Often correlates with nFeature_RNA."
      ),
      tags$li(
        strong("percent.mt:"), 
        " Percentage of mitochondrial reads. High values may indicate stressed or dying cells, so set an upper threshold where the majority of cells lie below."
      )
    ),
    hr(),
    
    # Guidance for filtering
    h3("How to Set Filters"),
    p("Use the violin plots for each metric to guide filtering:"),
    tags$ul(
      tags$li("Look for the bulk of cells (the dense 'mass' in the violin)."),
      tags$li("Exclude extreme outliers on either side."),
      tags$li("Aim to retain the majority of cells while removing low-quality ones.")
    ),
    p("After adjusting thresholds, click ", strong("Apply Filter"), " for that sample. The filtered dataset will automatically update."),
    hr(),
    
    # Step 1
    h3("Step 1: Inspect QC per Sample"),
    p("Adjust thresholds per sample below. Each panel includes its own Apply button."),
    uiOutput(ns("per_sample_panels")),
    
    hr(),
    # Step 2
    h3("Step 2: Review Summary"),
    DT::DTOutput(ns("qc_summary_table")),
    br(),
    verbatimTextOutput(ns("filter_info"))
  )
}


qcFilterServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    # ---------- helpers ----------
    `%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a
    .safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
    
    .choose_group_col <- function(obj) {
      if ("orig.ident" %in% colnames(obj@meta.data)) return("orig.ident")
      if ("sample" %in% colnames(obj@meta.data))     return("sample")
      stop("No 'orig.ident' or 'sample' column found in metadata.")
    }
    .guess_mito_pattern <- function(obj) {
      feats <- rownames(obj)
      if (length(feats) && any(grepl("^MT-", feats))) "^MT-" else "^mt-"
    }
    
    # Ensure percent.mt exists on the pooled object
    observeEvent(app_data$seurat_obj, {
      req(app_data$seurat_obj)
      obj <- app_data$seurat_obj
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        patt <- .guess_mito_pattern(obj)
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = patt)
        app_data$seurat_obj <- obj
      }
    }, ignoreInit = FALSE)
    
    # ---------- dynamic per-sample panels ----------
    output$per_sample_panels <- renderUI({
      req(app_data$seurat_obj)
      obj  <- app_data$seurat_obj
      grp  <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = grp)
      if (!length(parts)) return(tags$em("No samples found."))
      
      tagList(lapply(names(parts), function(s) {
        sid <- .safe_id(s)
        so  <- parts[[s]]
        
        # add percent.mt if missing in this split
        if (!"percent.mt" %in% colnames(so@meta.data)) {
          patt <- .guess_mito_pattern(so)
          so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = patt)
        }
        
        # default thresholds from quantiles
        qf_lo <- suppressWarnings(stats::quantile(so$nFeature_RNA, 0.05, na.rm = TRUE)) %||% 200
        qf_hi <- suppressWarnings(stats::quantile(so$nFeature_RNA, 0.95, na.rm = TRUE)) %||% 2500
        qmt   <- suppressWarnings(stats::quantile(so$percent.mt,  0.95, na.rm = TRUE)) %||% 5
        
        # per-sample violin (taller)
        local({
          .sid <- sid; .so <- so; .s <- s
          output[[paste0("vln_", .sid)]] <- renderPlot({
            VlnPlot(.so, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) +
              ggtitle(paste("QC violins —", .s))
          }, height = 600)
        })
        
        wellPanel(
          h4(sprintf("Sample: %s", s)),
          plotOutput(session$ns(paste0("vln_", sid)), height = "600px"),
          fluidRow(
            column(4, numericInput(session$ns(paste0("minF_",  sid)), "Min #features", round(qf_lo), min = 0)),
            column(4, numericInput(session$ns(paste0("maxF_",  sid)), "Max #features", round(qf_hi), min = 0)),
            column(4, numericInput(session$ns(paste0("maxMT_", sid)), "Max % mito", round(qmt, 1), min = 0, max = 100, step = 0.1))
          ),
          actionButton(session$ns(paste0("apply_", sid)), "Apply Filter for This Sample", class = "btn-primary")
        )
      }))
    })
    
    # ---------- attach one observer per Apply button (only once) ----------
    registered <- reactiveVal(character())
    
    observe({
      req(app_data$seurat_obj)
      obj  <- app_data$seurat_obj
      grp  <- .choose_group_col(obj)
      parts <- SplitObject(obj, split.by = grp)
      if (!length(parts)) return()
      
      sids <- sapply(names(parts), .safe_id)
      todo <- setdiff(sids, registered())
      if (!length(todo)) return()
      
      # map sid -> original sample name
      sid2name <- setNames(names(parts), sids)
      
      for (sid in todo) {
        btn_id <- paste0("apply_", sid)
        
        local({
          .sid <- sid
          observeEvent(input[[btn_id]], {
            req(app_data$seurat_obj)
            obj_now <- app_data$seurat_obj
            grp_now <- .choose_group_col(obj_now)
            parts_now <- SplitObject(obj_now, split.by = grp_now)
            
            s_name <- sid2name[[.sid]]
            if (is.null(s_name) || !(s_name %in% names(parts_now))) {
              showNotification("❌ Could not resolve sample for this button.", type = "error")
              return()
            }
            
            # read thresholds
            minF  <- input[[paste0("minF_",  .sid)]] %||% 200
            maxF  <- input[[paste0("maxF_",  .sid)]] %||% 2500
            maxMT <- input[[paste0("maxMT_", .sid)]] %||% 5
            
            showNotification(sprintf("🔧 Filtering sample '%s'…", s_name),
                             type = "message", id = "qc_progress", duration = NULL)
            
            so <- parts_now[[s_name]]
            if (!"percent.mt" %in% colnames(so@meta.data)) {
              patt <- .guess_mito_pattern(so)
              so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = patt)
            }
            
            filtered <- tryCatch({
              subset(
                so,
                subset = nFeature_RNA > minF &
                  nFeature_RNA < maxF &
                  percent.mt    < maxMT
              )
            }, error = function(e) {
              removeNotification("qc_progress")
              showNotification(paste("❌ Filtering error:", e$message), type = "error")
              return(NULL)
            })
            if (is.null(filtered)) return()
            
            # replace only this sample, keep others unchanged
            parts_now[[s_name]] <- filtered
            
            # drop empty samples (if any)
            nonempty <- sapply(parts_now, function(x) as.integer(ncol(x))) > 0
            if (!any(nonempty)) {
              removeNotification("qc_progress")
              showNotification("❌ All cells removed across samples. Loosen thresholds.", type = "error")
              return()
            }
            parts_now <- parts_now[nonempty]
            kept_names <- names(parts_now)
            
            # merge back (prefix barcodes by sample)
            merged <- if (length(parts_now) == 1) {
              one <- parts_now[[1]]
              colnames(one) <- paste0(kept_names[1], "_", colnames(one))
              one
            } else {
              merge(
                x            = parts_now[[1]],
                y            = parts_now[-1],
                add.cell.ids = kept_names,
                project      = "QC_filtered"
              )
            }
            
            app_data$seurat_obj <- merged
            
            # summary table
            cur_parts <- SplitObject(merged, split.by = "orig.ident")
            summary_df <- data.frame(
              sample = names(cur_parts),
              cells  = sapply(cur_parts, function(x) as.integer(ncol(x))),
              stringsAsFactors = FALSE
            )
            output$qc_summary_table <- DT::renderDT({
              DT::datatable(summary_df, options = list(pageLength = 10), rownames = FALSE)
            })
            output$filter_info <- renderText({
              paste("Total cells after filtering:", ncol(merged))
            })
            
            removeNotification("qc_progress")
            showNotification(sprintf("✅ Filter applied for sample: %s", s_name), type = "message")
          }, ignoreInit = TRUE)
        })
      }
      
      registered(c(registered(), todo))
    })
    
    # Initialize summary outputs
    output$qc_summary_table <- DT::renderDT({ NULL })
    output$filter_info <- renderText({ "" })
  })
}