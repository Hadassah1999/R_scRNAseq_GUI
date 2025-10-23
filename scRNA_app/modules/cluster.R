# tiny helper
`%||%` <- function(a, b) if (is.null(a)) b else a

clusterUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        sliderInput(
          ns("resolution"),
          "Resolution (higher = more clusters):",
          min = 0.1,
          max = 2.0,
          value = 0.5,
          step = 0.1
        ),
        actionButton(ns("run_clustering"), "Run Clustering"),
        br(), br(),
        helpText("Clusters are computed on PCA-reduced data."),
        br(),
        downloadButton(ns("download_plot"), "Download UMAP Plot"),
        br(), br(),
        helpText("After clustering, proceed to next step.")
      ),
      mainPanel(
        h4("UMAP Colored by Clusters"),
        withSpinner(plotOutput(ns("umap_plot")), type = 6),
        br(),
        # Cluster info (only renders after clustering_done() is TRUE)
        h4("Cluster Summary: Number of Cells per Cluster"),
        verbatimTextOutput(ns("cluster_info")),
        p("This table shows the number of cells in each cluster. Larger clusters usually represent major cell populations, whereas small clusters may be rare cell types or artifacts. Compare with the UMAP plot to interpret cluster structure.")
      )
    )
  )
}

clusterServer <- function(id, app_data, go_to_step5 = NULL) {
  moduleServer(id, function(input, output, session) {
    clustering_done <- reactiveVal(FALSE)
    
    observeEvent(input$run_clustering, {
      obj     <- app_data$seurat_obj %||% app_data$dataset
      num_pcs <- app_data$num_pcs
      
      if (is.null(obj) || is.null(num_pcs)) {
        showNotification("❌ Cannot run clustering. Make sure PCA (Step 3) was completed.", type = "error")
        return()
      }
      
      tryCatch({
        showNotification("ℹ️ Running clustering... Please wait", type = "message",
                         duration = NULL, id = "clustering_progress")
        
        # neighbors, clusters, UMAP
        obj <- FindNeighbors(obj, dims = 1:num_pcs)
        obj <- FindClusters(obj, resolution = input$resolution)
        Idents(obj) <- obj$seurat_clusters
        obj <- RunUMAP(obj, dims = 1:num_pcs)
        
        # persist useful params as bulletproof fallbacks
        if (is.null(obj@misc$pca))        obj@misc$pca        <- list()
        if (is.null(obj@misc$clustering)) obj@misc$clustering <- list()
        obj@misc$pca$pcs_used_for_neighbors <- num_pcs
        obj@misc$clustering$resolution       <- input$resolution
        
        # keep both references in sync so other modules (and summary) see updates
        app_data$seurat_obj <- obj
        app_data$dataset    <- obj
        
        clustering_done(TRUE)
        removeNotification("clustering_progress")
        showNotification("✅ Clustering completed!", type = "message")
        
        # optional: navigate to step 5 if your wizard uses it
        if (!is.null(go_to_step5)) go_to_step5()
        
      }, error = function(e) {
        removeNotification("clustering_progress")
        showNotification(paste("❌ Clustering error:", e$message), type = "error")
      })
    })
    
    # UMAP plot with cluster labels
    output$umap_plot <- renderPlot({
      req(clustering_done())
      obj <- app_data$seurat_obj %||% app_data$dataset
      req(obj)
      req("umap" %in% names(obj@reductions))
      
      DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
        ggtitle("UMAP: Clusters")
    })
    
    # Cluster info table
    output$cluster_info <- renderPrint({
      req(clustering_done())
      obj <- app_data$seurat_obj %||% app_data$dataset
      req(obj)
      cat("Number of cells per cluster:\n")
      print(table(Idents(obj)))
    })
    
    # Download UMAP plot
    output$download_plot <- downloadHandler(
      filename = function() paste0("UMAP_Clusters_", Sys.Date(), ".png"),
      content = function(file) {
        obj <- app_data$seurat_obj %||% app_data$dataset
        req(obj)
        ggsave(
          file,
          plot = DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
            ggtitle("UMAP: Clusters"),
          width = 7, height = 5
        )
      }
    )
    
    # RESET functionality (kept, but now clears both references)
    observeEvent(app_data$reset_cluster, {
      updateSliderInput(session, "resolution", value = 0.5)
      if (!is.null(app_data$seurat_obj)) {
        obj <- app_data$seurat_obj
        if ("seurat_clusters" %in% colnames(obj@meta.data)) obj$seurat_clusters <- NULL
        obj@reductions$umap <- NULL
        obj@graphs <- list()
        app_data$seurat_obj <- obj
        app_data$dataset    <- obj
      }
      clustering_done(FALSE)
      app_data$reset_cluster <- FALSE
    })
  })
}