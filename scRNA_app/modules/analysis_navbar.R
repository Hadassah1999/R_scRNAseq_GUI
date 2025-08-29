source("modules/cellAnnotation.R")
source("modules/gene_highlighting.R")
source("modules/data_summary.R")
source("modules/marker_genes.R")
source("modules/download.R")
analysisUI <- function(id) {
  ns <- NS(id)
  navbarPage(
    title = "Analysis",
    id = ns("navbar"),
    
    header = tags$div(
      style = "position: absolute; right: 15px; top: 8px; z-index: 1050;",
      actionButton(ns("back_btn"), "⬅ Back to Upload")
    ),
    tabPanel("Summary", summaryUI(ns("summary"))),
    tabPanel("Cell Annotation", annotationUI(ns("annotation"))),
    tabPanel("Cell highlighting by gene", geneHighlightingUI(ns("geneHighlighting"))),
    tabPanel("Marker Genes", markersUI(ns("markers"))),
    tabPanel("Download rds", downloadUI(ns("download")))

  )
}

analysisServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    summaryServer("summary", app_data)
    annotationServer("annotation", app_data)
    geneHighlightingServer("geneHighlighting", app_data)
    markersServer("markers", app_data)
    downloadServer("download", app_data)
    observeEvent(input$back_btn, {
      session$sendCustomMessage(type = "clear_notifications", message = list())
      app_data$file_type <- NULL
      app_data$dataset <- NULL
      app_data$proceed <- FALSE
      app_data$reset_upload <- TRUE
      app_data$reset_qc <- TRUE
      app_data$reset_cluster <- TRUE
      app_data$reset_summary <- TRUE
      app_data$reset_final <- TRUE
      app_data$reset_gene_highlighting <- TRUE
      app_data$reset_marker_module <- TRUE
      app_data$reset_normalization <- TRUE
      app_data$reset_pca <- TRUE
      app_data$reset_filter <- TRUE
    })
  })
}