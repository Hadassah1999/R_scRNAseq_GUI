source("modules/qc.R")
source("modules/visualize.R")
source("modules/normalize.R")
source("modules/pca.R")
source("modules/cluster.R")
source("modules/cellAnnotation.R")
source("modules/cellEval.R")
source("modules/gene_highlighting.R")
source("modules/marker_genes.R")
source("modules/final.R")
source("modules/pathways.R")
source("modules/cell_subtypes.R")



qcUI <- function(id) {
  ns <- NS(id)
  navbarPage(
    "Preprocessing / QC",
    id = ns("qc_navbar"),
    
    tabPanel("Present data", qcModuleUI(ns("qc_step"))),
    header = tags$div(
      style = "position: absolute; right: 15px; top: 8px; z-index: 1050;",
      actionButton(ns("back_btn"), "⬅ Back to Upload")
    ),
    
    tabPanel("Filter data", 
             qcFilterUI(ns("filter_step"))
    ),
    
    tabPanel("Normalization", 
             normalizationUI(ns("norm_step"))
    ),
    
    tabPanel("PCA", pcaUI(ns("pca_step"))),
    
    tabPanel("Clustering", clusterUI(ns("cluster_step"))),
    
    tabPanel("Cell Annotation", annotationUI(ns("ann_step"))),
    
    tabPanel("Cell Evaluation", evalUI(ns("Evaluation"))),
    
    tabPanel("Cell highlighting by gene", geneHighlightingUI(ns("geneHighlighting"))),
    
    tabPanel("Pathway Analysis", pathwayUI(ns("pathway"))),
    
    tabPanel("Cell Subtypes", subtypeUI(ns("subtypes"))),
    
    tabPanel("Marker Genes", markersUI(ns("markers"))),
    
    tabPanel("Export", finalUI(ns("final_step")))
  )
}

qcServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    qcModuleServer("qc_step", app_data)
    qcFilterServer("filter_step", app_data)
    normalizationServer("norm_step", app_data)
    pcaServer("pca_step", app_data)
    clusterServer("cluster_step", app_data)
    annotationServer("ann_step", app_data)
    evalServer("Evaluation", app_data)
    geneHighlightingServer("geneHighlighting", app_data)
    markersServer("markers", app_data)
    subtypeServer("subtypes", app_data)
    pathwayServer("pathway", app_data)
    finalServer("final_step", app_data)
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