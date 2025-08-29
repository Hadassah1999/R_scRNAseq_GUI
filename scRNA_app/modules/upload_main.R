uploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(
      tags$style(HTML("
        .big-radio .shiny-options-group label { font-size: 18px; }
      "))
    ),
    div(
      style = "font-size: 45px; font-weight: bold; margin-top: 30px; margin-bottom: 20px;",
      "🔬 Upload Your scRNA-seq Dataset: "
    ),
    fluidRow(
      column(
        width = 6, offset = 3,
        style = "margin-top: 40px;",
        wellPanel(
          h3("Step 1: Choose Data Format"),
          div(
            class = "big-radio",
            radioButtons(
              inputId = ns("file_type"),
              label = tags$span("Select the type of input data:", style = "font-size: 20px;"),
              choices = c("Seurat (.rds)" = "rds", "10X Directory" = "10x", "HDF5 (.h5)" = "h5"),
              selected = "rds"
            )
          ),
          
          hr(),
          
          # Step 2 is dynamic: show path only for RDS, otherwise show an info note
          uiOutput(ns("path_block_ui")),
          
          hr(),
          
          actionButton(ns("submit_btn"), "📂 Load Data", class = "btn-primary", style = "font-size: 18px;"),
          br(), br(),
          
          textOutput(ns("status")),
          br(),
          
          uiOutput(ns("continue_btn_ui"))
        )
      )
    )
  )
}




uploadServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    # Keep dataset reactive here so we can preview summary (if you use it)
    dataset <- reactiveVal(app_data$dataset)
    
    # Dynamic Step 2 block
    output$path_block_ui <- renderUI({
      if (identical(input$file_type, "rds")) {
        tagList(
          h3("Step 2: Provide File Path"),
          textInput(
            inputId = session$ns("file_path"),
            label = tags$span("Enter the full path to your dataset:", style = "font-size: 20px;"),
            placeholder = "/Users/username/Documents/my_data.rds"
          )
        )
      } else {
        tagList(
          h3("Step 2: No Path Needed"),
          p("For 10X/HDF5 inputs, you'll add sample paths in the QC step.")
        )
      }
    })
    
    # Clear stale UI when changing file type
    observeEvent(input$file_type, {
      output$status <- renderText("")
      output$continue_btn_ui <- renderUI(NULL)
      if (identical(input$file_type, "rds")) {
        updateTextInput(session, "file_path", value = "")
      }
    }, ignoreInit = TRUE)
    
    # Optional: summary renderer if you show it somewhere
    output$seurat_summary <- renderText({
      req(dataset())
      obj <- dataset()
      assays <- tryCatch(paste(Assays(obj), collapse = ", "), error = function(e) "Unknown")
      reductions <- tryCatch(paste(names(obj@reductions), collapse = ", "), error = function(e) "None")
      identities <- tryCatch(length(unique(Idents(obj))), error = function(e) NA)
      default_assay <- tryCatch(DefaultAssay(obj), error = function(e) "Unknown")
      clusters_detected <- if (!is.na(identities)) paste0(identities, " cluster(s)") else "Not clustered"
      paste0(
        "🧬 Class: ", class(obj)[1], "\n",
        "🔢 Cells: ", ncol(obj), "\n",
        "📈 Features: ", nrow(obj), "\n",
        "🧪 Default Assay: ", default_assay, "\n",
        "🧩 Assays: ", assays, "\n",
        "📉 Reductions: ", reductions, "\n",
        "🧬 Clusters: ", clusters_detected
      )
    })
    
    observeEvent(input$submit_btn, {
      ft <- input$file_type
      
      # RDS path flow
      if (identical(ft, "rds")) {
        path <- input$file_path
        if (!nzchar(path) || !file.exists(path)) {
          output$status <- renderText("❌ File or directory does not exist at the given path.")
          output$continue_btn_ui <- renderUI(NULL)
          return()
        }
        
        app_data$file_type <- "rds"
        tryCatch({
          obj <- readRDS(path)
          if (!inherits(obj, "Seurat")) {
            showNotification("❌ The file does not contain a Seurat object.", type = "error")
            output$status <- renderText("❌ Not a Seurat object.")
            output$continue_btn_ui <- renderUI(NULL)
            return()
          }
          dataset(obj)
          app_data$dataset <- obj
          output$status <- renderText("✅ Seurat object successfully loaded.")
          showNotification("✅ File loaded successfully!", type = "message")
          output$continue_btn_ui <- renderUI({
            actionButton(session$ns("continue_btn"), "➡ Continue to Analysis", class = "btn-success")
          })
        }, error = function(e) {
          output$status <- renderText(paste("❌ Error loading file:", e$message))
          output$continue_btn_ui <- renderUI(NULL)
          showNotification("❌ Error: Check your file path and format.", type = "error")
        })
        
      } else {
        # 10X / H5 flow — no path required here
        app_data$file_type <- ft
        app_data$raw_file_path <- NULL  # ensure no stale path
        output$status <- renderText("✅ Format selected. You'll add sample paths in the QC step.")
        showNotification("✅ Format accepted! Proceed to preprocessing.", type = "message")
        output$continue_btn_ui <- renderUI({
          actionButton(session$ns("continue_btn"), "➡ Continue to Analysis", class = "btn-success")
        })
      }
    })
    
    observeEvent(input$continue_btn, {
      showNotification("📊 Proceeding to next step...", type = "message")
      app_data$proceed <- TRUE
    })
    
    # RESET hook
    observeEvent(app_data$reset_upload, {
      dataset(NULL)
      app_data$dataset <- NULL
      app_data$file_type <- NULL
      app_data$raw_file_path <- NULL
      app_data$proceed <- FALSE
      updateRadioButtons(session, "file_type", selected = "rds")
      # file_path input is inside path_block_ui; it will be cleared on re-render
      output$status <- renderText("")
      output$continue_btn_ui <- renderUI(NULL)
      app_data$reset_upload <- FALSE
    })
  })
}
