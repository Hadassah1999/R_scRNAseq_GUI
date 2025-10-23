library(shiny)

# Source global settings and modules
source("global.R")
source("modules/upload_main.R")
source("modules/analysis_navbar.R")
source("modules/qc_navbar.R")

options(shiny.legacy.datatable = FALSE)

# Store user data and app state
app_data <- reactiveValues(
  file_type = NULL,
  dataset = NULL,
  proceed = FALSE
)

# Define UI
ui <- fluidPage(
  tags$head(
    tags$title("scRNAseq App"),
    tags$link(rel = "icon", type = "image/png", href = "favicon6.png"),
    tags$script(HTML("
    Shiny.addCustomMessageHandler('clear_notifications', function(message) {
      $('.shiny-notification').remove();
    });
  "))
  ),
  uiOutput("dynamic_ui")
)

app_data <- reactiveValues(
  file_type = NULL,
  dataset = NULL,
  proceed = FALSE
)

# Define Server
server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # Track which module has been loaded
  loaded_modules <- reactiveValues(
    upload = FALSE,
    analysis = FALSE,
    qc = FALSE
  )
  
  # Dynamically load only the needed module server
  observe({
    if (!app_data$proceed && !loaded_modules$upload) {
      uploadServer("upload", app_data)
      loaded_modules$upload <- TRUE
    }
    
    if (app_data$proceed && app_data$file_type == "rds" && !loaded_modules$analysis) {
      analysisServer("analysis", app_data)
      loaded_modules$analysis <- TRUE
    }
    
    if (app_data$proceed && app_data$file_type %in% c("10x", "h5") && !loaded_modules$qc) {
      qcServer("qc", app_data)
      loaded_modules$qc <- TRUE
    }
  })
  
  # Render UI based on state
  output$dynamic_ui <- renderUI({
    if (!app_data$proceed) {
      uploadUI("upload")
    } else if (app_data$file_type == "rds") {
      analysisUI("analysis")
    } else if (app_data$file_type %in% c("10x", "h5")) {
      qcUI("qc")
    }
  })
}




# Launch the app
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
