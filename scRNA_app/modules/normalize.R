normalizationUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("norm_method"), "Normalization method:",
                choices = c("LogNormalize", "CLR", "RC"),
                selected = "LogNormalize"),
    numericInput(ns("scale_factor"), "Scale factor (only for LogNormalize):", value = 10000, min = 1),
    checkboxInput(ns("do_scale"), "Also scale data", value = TRUE),
    actionButton(ns("apply_norm"), "Apply Normalization"),
    verbatimTextOutput(ns("norm_status"))
  )
}

normalizationServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$apply_norm, {
      req(app_data$seurat_obj)
      tryCatch({
        app_data$seurat_obj <- NormalizeData(
          app_data$seurat_obj,
          normalization.method = input$norm_method,
          scale.factor = if (input$norm_method == "LogNormalize") input$scale_factor else NULL
        )
        app_data$seurat_obj <- FindVariableFeatures(app_data$seurat_obj)
        if (input$do_scale) app_data$seurat_obj <- ScaleData(app_data$seurat_obj)
        output$norm_status <- renderText({
          paste0("✅ Normalization done using method: ", input$norm_method,
                 ifelse(input$do_scale, " and scaling applied.", " (scaling skipped)."),
                 "\n\nProceed to next step.")
        })
      }, error = function(e) {
        output$norm_status <- renderText(paste("❌ Error:", e$message))
      })
    })
    
    # RESET
    observeEvent(app_data$reset_normalization, {
      updateSelectInput(session, "norm_method", selected = "LogNormalize")
      updateNumericInput(session, "scale_factor", value = 10000)
      updateCheckboxInput(session, "do_scale", value = TRUE)
      output$norm_status <- renderText("")   # safe re-render of plain text
      app_data$reset_normalization <- FALSE
    })
  })
}

