
downloadUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Save Your Results"),
    p("Download the current Seurat object (.rds) to preserve your work for future analysis."),
    textInput(ns("filename"), "File name (without extension):", value = "processed_data"),
    downloadButton(ns("download_rds"), "Download .rds", class = "btn-primary"),
    br(), br(),
    textOutput(ns("status"))
  )
}

downloadServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0 || identical(a, "")) b else a
    
    # Disable button initially
    observe({
      shinyjs::disable("download_rds")
    })
    
    # Enable button when object is available
    observe({
      if (!is.null(app_data$dataset)) {
        shinyjs::enable("download_rds")
      } else {
        shinyjs::disable("download_rds")
      }
    })
    
    output$download_rds <- downloadHandler(
      filename = function() {
        paste0(input$filename %||% "processed_data", ".rds")
      },
      content = function(file) {
        if (is.null(app_data$dataset)) {
          output$status <- renderText("⚠️ No Seurat object available to save.")
          return()
        }
        tryCatch({
          saveRDS(app_data$dataset, file)
          output$status <- renderText(paste0("✅ File saved as ", input$filename, ".rds"))
        }, error = function(e) {
          output$status <- renderText(paste0("❌ Error saving file: ", e$message))
        })
      }
    )
    
  })
}
