pcaUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("npc"), "Number of PCs to compute:", value = 30, min = 5, max = 50),
        actionButton(ns("run_pca"), "Run PCA"),
        br(),
        helpText("Step 1: Choose how many principal components (PCs) to compute."),
        helpText("Step 2: Use the Elbow Plot to decide how many PCs capture most variance."),
        br(),
        uiOutput(ns("pc_selector_ui")),
        br(),
        helpText("You can now proceed to step 4.")
      ),
      mainPanel(
        h4("Elbow Plot"),
        withSpinner(plotOutput(ns("elbow_plot")), type = 6),
        h4("PCA Scatter Plot"),
        withSpinner(plotOutput(ns("pca_scatter")), type = 6)
      )
    )
  )
}


pcaServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    pca_result  <- reactiveVal(NULL)
    pcs_selected<- reactiveVal(NULL)
    
    observeEvent(input$run_pca, {
      obj <- app_data$seurat_obj
      if (is.null(obj)) {
        showNotification("❌ No Seurat object found. Please load and preprocess your data first.", type="error")
        return()
      }
      tryCatch({
        showNotification("ℹ️ Running PCA... Please wait.", type="message", duration=NULL, id="pca_running")
        obj <- RunPCA(obj, npcs = input$npc)
        removeNotification("pca_running")
        showNotification("✅ PCA completed.", type="message")
        app_data$seurat_obj <- obj
        pca_result(obj)
        pcs_selected(10)
      }, error = function(e) {
        removeNotification("pca_running")
        showNotification(paste("❌ PCA error:", e$message), type="error")
      })
    })
    
    output$elbow_plot <- renderPlot({
      req(pca_result())
      ElbowPlot(pca_result(), ndims = input$npc)
    })
    
    output$pca_scatter <- renderPlot({
      req(pca_result())
      DimPlot(pca_result(), reduction="pca", dims=c(1,2)) + ggtitle("PCA Scatter Plot (PC1 vs PC2)")
    })
    
    output$pc_selector_ui <- renderUI({
      req(pca_result())
      numericInput(session$ns("pcs_to_use"), "Number of PCs to use in next step:",
                   value = pcs_selected() %||% 10, min = 2, max = input$npc)
    })
    
    observeEvent(input$pcs_to_use, {
      pcs_selected(input$pcs_to_use)
      app_data$num_pcs <- input$pcs_to_use
    })
    
    # RESET
    observeEvent(app_data$reset_pca, {
      pca_result(NULL)
      pcs_selected(NULL)
      updateNumericInput(session, "npc", value = 30)
      # optional: clear pca reduction from object
      if (!is.null(app_data$seurat_obj)) {
        obj <- app_data$seurat_obj
        obj@reductions$pca <- NULL
        app_data$seurat_obj <- obj
      }
      app_data$reset_pca <- FALSE
    })
  })
}
