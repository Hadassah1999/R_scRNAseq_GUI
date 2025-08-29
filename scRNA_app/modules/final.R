finalUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Processing Complete!"),
    p("Your .rds file has been created. Name it and download it below."),
    textInput(ns("filename"), "Name your .rds file (without extension):", value = "processed_data"),
    downloadButton(ns("download_rds"), "Download Processed .rds"),
    br(), br(),
    # Placeholder for the instruction message (shown dynamically)
    uiOutput(ns("next_step_msg"))
  )
}

finalServer <- function(id, app_data, proceed_callback = NULL) {
  moduleServer(id, function(input, output, session) {
    `%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0 || identical(a, "")) b else a
    
    rds_downloaded <- reactiveVal(FALSE)
    
    # Handle file download
    output$download_rds <- downloadHandler(
      filename = function() paste0(input$filename %||% "processed_data", ".rds"),
      content  = function(file) {
        if (is.null(app_data$seurat_obj)) stop("No Seurat object found to save.")
        saveRDS(app_data$seurat_obj, file)
        rds_downloaded(TRUE)
      }
    )
    
    # Show message after successful download
    output$next_step_msg <- renderUI({
      req(rds_downloaded())
      tags$p(
        style = "font-size: 16px; color: #444; margin-top: 20px;",
        "To proceed to the analysis step, return to the Upload page and upload the ",
        strong("Seurat (.rds) file"), " you just created."
      )
    })
    
    # Reset hook
    observeEvent(app_data$reset_final, {
      updateTextInput(session, "filename", value = "processed_data")
      rds_downloaded(FALSE)
      output$next_step_msg <- renderUI(NULL)
      app_data$reset_final <- FALSE
    })
  })
}