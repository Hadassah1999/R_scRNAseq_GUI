summaryUI <- function(id) {
  ns <- NS(id)
  tagList(
    mainPanel(
      h4("Dataset Summary"),
      actionButton(ns("refresh"), "Refresh Summary"),
      br(),
      verbatimTextOutput(ns("seurat_summary")),
      br(),
      h4("Metadata Preview"),
      tableOutput(ns("meta_table"))
    )
  )
}
.extract_clust_params <- function(obj) {
  pcs_txt <- "Unknown"
  res_txt <- "Unknown"
  
  cmds <- tryCatch(obj@commands, error = function(e) NULL)
  if (!is.null(cmds) && length(cmds) > 0) {
    # names can be "FindNeighbors", "FindNeighbors.RNA", "FindNeighbors.SCT", etc.
    fn_names <- grep("^FindNeighbors", names(cmds), value = TRUE, perl = TRUE)
    fc_names <- grep("^FindClusters", names(cmds), value = TRUE, perl = TRUE)
    
    # --- PCs from FindNeighbors params ---
    if (length(fn_names) > 0) {
      fn <- cmds[[tail(fn_names, 1)]]
      dims <- NULL
      # Seurat v4/v5: prefer @params; fallback to @call
      dims <- tryCatch(fn@params$dims, error = function(e) NULL)
      if (is.null(dims)) dims <- tryCatch(fn@call$dims, error = function(e) NULL)
      if (is.numeric(dims)) pcs_txt <- as.character(max(dims))
    }
    
    # --- Resolution from FindClusters params ---
    if (length(fc_names) > 0) {
      fc <- cmds[[tail(fc_names, 1)]]
      res <- tryCatch(fc@params$resolution, error = function(e) NULL)
      if (is.null(res)) res <- tryCatch(fc@call$resolution, error = function(e) NULL)
      if (!is.null(res)) res_txt <- as.character(res)
    }
  }
  
  # Fallbacks from @misc if command log is missing/stripped
  if (identical(pcs_txt, "Unknown")) {
    pcs_txt <- tryCatch(as.character(obj@misc$pca$pcs_used_for_neighbors), error = function(e) "Unknown")
  }
  if (identical(res_txt, "Unknown")) {
    res_txt <- tryCatch(as.character(obj@misc$clustering$resolution), error = function(e) "Unknown")
  }
  
  list(pcs = pcs_txt, resolution = res_txt)
}

# ----------------------------
# SUMMARY SERVER (updated)
# ----------------------------
summaryServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    refresh_tick <- reactiveVal(0)
    observeEvent(input$refresh, { refresh_tick(refresh_tick() + 1) })
    
    output$seurat_summary <- renderText({
      refresh_tick()
      obj <- app_data$seurat_obj %||% app_data$dataset
      req(obj)
      
      assays         <- tryCatch(paste(Assays(obj), collapse = ", "), error = function(e) "Unknown")
      reductions     <- tryCatch(paste(names(obj@reductions), collapse = ", "), error = function(e) "None")
      identities_n   <- tryCatch(length(unique(Idents(obj))), error = function(e) NA)
      default_assay  <- tryCatch(DefaultAssay(obj), error = function(e) "Unknown")
      clusters_desc  <- if (!is.na(identities_n)) paste0(identities_n, " cluster(s)") else "Not clustered"
      
      cp <- .extract_clust_params(obj)
      
      paste0(
        "📦 Class: ", class(obj)[1], "\n",
        "🧬 Default Assay: ", default_assay, "\n",
        "🧪 Assays: ", assays, "\n",
        "📉 Reductions: ", reductions, "\n",
        "🧫 Cells: ", ncol(obj), "\n",
        "🧬 Features: ", nrow(obj), "\n",
        "🗂 Metadata Columns: ", ncol(obj@meta.data), "\n",
        "🔍 Clustering: ", clusters_desc, "\n",
        "⚙️ PCs: ", cp$pcs, "\n",
        "⚙️ Resolution: ", cp$resolution
      )
    })
    outputOptions(output, "seurat_summary", suspendWhenHidden = FALSE)
    
    output$meta_table <- renderTable({
      refresh_tick()
      obj <- app_data$seurat_obj %||% app_data$dataset
      req(obj)
      head(obj@meta.data, 10)
    })
    outputOptions(output, "meta_table", suspendWhenHidden = FALSE)
    
    observeEvent(app_data$reset_summary, {
      app_data$reset_summary <- FALSE
    })
  })
}