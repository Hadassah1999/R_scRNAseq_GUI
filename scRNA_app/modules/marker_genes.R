
markersUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Clustering & Marker Finding",
    sidebarLayout(
      sidebarPanel(
        helpText(
          "This page identifies marker genes for each cluster and shows the top markers you choose.",
          "⚠️ The object must already be clustered (i.e., have a column with cluster identities).",
          "By default, the column 'seurat_clusters' is used. Change below if your clustering is stored under another name."
        ),
        tags$hr(),
        textInput(ns("cluster_col"), "Cluster column name in metadata:", value = "seurat_clusters"),
        numericInput(ns("top_n"), "Top markers per cluster:", value = 2, min = 1),
        actionButton(ns("run_all_btn"), "Find Markers"),
        br(), br(),
        downloadButton(ns("download_markers"), "Download Top-Markers Table"),
        br(), br(),
        helpText("Note: On large datasets this step can take several minutes.")
      ),
      mainPanel(
        tags$style(HTML(
          ".center-wrap{display:flex;align-items:center;justify-content:center;min-height:320px;}
           .progress-bar-outer{width:60%;max-width:520px;background:#f3f3f3;border:1px solid #ddd;border-radius:6px;}
           .progress-bar-inner{width:100%;height:18px;background:#4CAF50;border-radius:6px;}
           .muted{color:#666;font-size:0.95em;}"
        )),
        uiOutput(ns("dynamic_panel"))
      )
    )
  )
}

markersServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    calculating <- reactiveVal(FALSE)
    `%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a
    
    .order_vec <- function(df) {
      if ("avg_log2FC" %in% names(df)) return(df$avg_log2FC)
      if ("avg_logFC"  %in% names(df)) return(df$avg_logFC)
      if ("p_val_adj"  %in% names(df)) return(-df$p_val_adj)
      if ("p_val"      %in% names(df)) return(-df$p_val)
      rep(0, nrow(df))
    }
    
    output$dynamic_panel <- renderUI({
      ns <- session$ns
      if (isTRUE(calculating())) {
        div(class = "center-wrap",
            div(
              div(h4("Processing…")),
              p(class = "muted",
                "Marker identification is running. This may take a few minutes depending on your dataset size."),
              br(),
              div(class = "progress-bar-outer",
                  div(class = "progress-bar-inner"))
            )
        )
      } else {
        fluidRow(
          column(
            width = 8,
            h4("Top Markers by Cluster"),
            tableOutput(ns("top_markers_preview"))
          ),
          column(
            width = 4,
            uiOutput(ns("col_explainer"))
          )
        )
      }
    })
    
    output$col_explainer <- renderUI({
      req(!isTRUE(calculating()))
      mk <- app_data$top_markers
      req(mk)
      cols <- colnames(mk)
      desc <- list(
        gene       = "Gene symbol (rownames promoted to a column).",
        cluster    = "Cluster ID the marker was computed for.",
        p_val      = "Raw p-value for differential expression.",
        p_val_adj  = "Adjusted p-value (multiple testing correction).",
        avg_log2FC = "Average log2 fold-change (cluster vs others). Higher is more upregulated.",
        avg_logFC  = "Average log fold-change (legacy column name in some Seurat versions).",
        pct.1      = "Fraction of cells in the cluster expressing the gene.",
        pct.2      = "Fraction of cells outside the cluster expressing the gene.",
        power      = "Test statistic / power (present in some methods)."
      )
      present <- intersect(names(desc), cols)
      if (!length(present)) return(NULL)
      tags$div(
        tags$hr(),
        tags$h5("Columns explained"),
        tags$ul(lapply(present, function(nm) tags$li(tags$b(nm), ": ", desc[[nm]])))
      )
    })
    
    observeEvent(input$run_all_btn, {
      obj <- app_data$seurat_obj %||% app_data$dataset
      req(obj)
      
      # --- get cluster column name from input ---
      clcol <- input$cluster_col %||% "seurat_clusters"
      
      # validate column exists and has at least two groups
      if (!(clcol %in% colnames(obj@meta.data)) ||
          length(unique(obj@meta.data[[clcol]])) < 2) {
        showNotification(sprintf("❌ Column '%s' not found or has <2 clusters.", clcol), type = "error")
        return()
      }
      
      if (is.na(input$top_n) || input$top_n < 1) {
        showNotification("Please provide a valid Top markers per cluster (>=1).", type = "error")
        return()
      }
      
      calculating(TRUE)
      showNotification("🔄 Finding markers for existing clusters…", type = "message",
                       duration = NULL, id = "progress")
      
      tryCatch({
        obj2 <- obj
        DefaultAssay(obj2) <- "RNA"
        if ("layers" %in% slotNames(obj2@assays$RNA)) {
          obj2 <- JoinLayers(obj2)
        }
        
        Idents(obj2) <- obj2@meta.data[[clcol]]
        
        markers <- FindAllMarkers(
          obj2,
          only.pos        = TRUE,
          min.pct         = 0.25,
          logfc.threshold = 0.25
        )
        app_data$all_markers <- markers
        
        if (!is.null(markers) && nrow(markers) > 0) {
          clusters <- sort(unique(as.character(markers$cluster)))
          top_list <- lapply(clusters, function(cl) {
            mk <- markers[markers$cluster == cl, , drop = FALSE]
            if (nrow(mk) == 0) return(NULL)
            mk <- mk[order(.order_vec(mk), decreasing = TRUE, na.last = NA), , drop = FALSE]
            mk[seq_len(min(input$top_n, nrow(mk))), , drop = FALSE]
          })
          app_data$top_markers <- dplyr::bind_rows(top_list[!vapply(top_list, is.null, logical(1))])
          rownames(app_data$top_markers) <- NULL
        } else {
          app_data$top_markers <- data.frame(message = "No markers were identified.")
        }
        
        calculating(FALSE)
        removeNotification("progress")
        showNotification(sprintf("✅ Markers computed using column '%s'.", clcol), type = "message")
        
      }, error = function(e) {
        calculating(FALSE)
        removeNotification("progress")
        showNotification(paste("❌ Error:", conditionMessage(e)), type = "error")
        app_data$top_markers <- data.frame(message = conditionMessage(e))
      })
    })
    
    observeEvent(input$top_n, {
      mk <- app_data$all_markers
      req(mk)
      if (nrow(mk) == 0) { app_data$top_markers <- mk; return() }
      clusters <- sort(unique(as.character(mk$cluster)))
      top_list <- lapply(clusters, function(cl) {
        mk_cl <- mk[mk$cluster == cl, , drop = FALSE]
        if (nrow(mk_cl) == 0) return(NULL)
        mk_cl <- mk_cl[order(.order_vec(mk_cl), decreasing = TRUE, na.last = NA), , drop = FALSE]
        mk_cl[seq_len(min(input$top_n, nrow(mk_cl))), , drop = FALSE]
      })
      top_tbl <- dplyr::bind_rows(top_list[!vapply(top_list, is.null, logical(1))])
      rownames(top_tbl) <- NULL
      app_data$top_markers <- top_tbl
    })
    
    output$top_markers_preview <- renderTable({
      req(app_data$top_markers)
      as.data.frame(app_data$top_markers)
    }, rownames = FALSE)
    outputOptions(output, "top_markers_preview", suspendWhenHidden = FALSE)
    
    output$download_markers <- downloadHandler(
      filename = function() paste0("top_markers_by_cluster_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(app_data$top_markers)
        utils::write.csv(as.data.frame(app_data$top_markers), file, row.names = FALSE)
      }
    )
    
    observeEvent(app_data$reset_marker_module, {
      app_data$all_markers <- NULL
      app_data$top_markers <- NULL
      calculating(FALSE)
      updateTextInput(session, "cluster_col", value = "seurat_clusters")
      updateNumericInput(session, "top_n", value = 2)
      app_data$reset_marker_module <- FALSE
    })
  })
}
