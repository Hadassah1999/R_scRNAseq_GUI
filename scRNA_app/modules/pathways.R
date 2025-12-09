# -----------------------------
# UI
# -----------------------------
pathwayUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarPanel(
      radioButtons(
        ns("species"),
        "Select species:",
        choices = c("Human", "Mouse"),
        selected = "Human"
      ),
      
      selectInput(
        ns("select_pathway"),
        "Select Pathway:",
        choices = c("cGAS-STING", "p53", "RIG-I/MDA5", "ISG", "Custom"),
        selected = "Custom"
      ),
      
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Custom'", ns("select_pathway")),
        textAreaInput(ns("custom_genes"), "Enter gene list (comma or newline separated):")
      ),
      
      actionButton(ns("apply_pathway"), "Run Pathway Scoring"),
      hr(),
      
      uiOutput(ns("celltype_selector_ui"))
    ),
    
    mainPanel(
      h4("Average pathway score by dataset"),
      plotOutput(ns("metadata_barplot")),
      tableOutput(ns("metadata_table"))
    )
  )
}


# -----------------------------
# SERVER
# -----------------------------
pathwayServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    pathway_genes_human <- list(
      "cGAS-STING" = c("CGAS", "TMEM173", "TBK1", "IRF3", "IRF7"),
      "p53"        = c("TP53", "CDKN1A", "GADD45A", "BAX", "MDM2"),
      "RIG-I/MDA5" = c("DDX58", "IFIH1", "IRF3", "IRF7", "ISG15"),
      "ISG"        = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2")
    )
    
    pathway_genes_mouse <- list(
      "cGAS-STING" = c("Mb21d1", "Tmem173", "Tbk1", "Irf3", "Irf7"),
      "p53"        = c("Trp53", "Cdkn1a", "Gadd45a", "Bax", "Mdm2"),
      "RIG-I/MDA5" = c("Ddx58", "Ifih1", "Irf3", "Irf7", "Isg15"),
      "ISG"        = c("Isg15", "Ifit1", "Ifit3", "Mx1", "Oas1a", "Oas2")
    )
    
    current_gene_set <- reactive({
      req(input$select_pathway)
      if (input$select_pathway %in% names(pathway_genes_human)) {
        if (input$species == "Mouse") pathway_genes_mouse[[input$select_pathway]]
        else                         pathway_genes_human[[input$select_pathway]]
      } else {
        genes <- unlist(strsplit(input$custom_genes, "[,\\n\\r]+"))
        trimws(genes[genes != ""])
      }
    })
    
    avg_scores_reactive <- reactiveVal(NULL)
    
    observeEvent(input$apply_pathway, {
      seurat_obj <- app_data$dataset %||% app_data$seurat_obj
      req(is(seurat_obj, "Seurat"))
      
      genes <- current_gene_set()
      genes_in_data <- genes[genes %in% rownames(seurat_obj)]
      
      if (length(genes_in_data) == 0) {
        showNotification("None of the selected genes are found in the dataset.", type = "error")
        return(NULL)
      }
      
      notif <- showNotification("Calculating UCell score…", duration = NULL, type = "message")
      on.exit(removeNotification(notif))
      
      # Run UCell
      gene_sets <- list(Pathway = genes_in_data)
      seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, gene_sets)
      score_col <- "Pathway_UCell"
      
      if (!score_col %in% colnames(seurat_obj@meta.data)) {
        showNotification("UCell scoring failed.", type="error")
        return(NULL)
      }
      
      if ("reactiveVal" %in% class(app_data)) app_data(seurat_obj)
      
      meta <- seurat_obj@meta.data
      if (!all(c("orig.ident", "cell_ann", score_col) %in% colnames(meta))) {
        showNotification("Missing orig.ident or cell_ann metadata.", type="error")
        return(NULL)
      }
      
      # -----------------------------
      # Average scores per origin × annotation + num cells + num positive
      # -----------------------------
      avg_scores <- meta %>%
        dplyr::select(orig.ident, cell_ann, !!score_col) %>%
        dplyr::rename(origin = orig.ident,
                      annotation = cell_ann,
                      score = !!score_col) %>%
        dplyr::group_by(origin, annotation) %>%
        dplyr::summarize(
          avg_score = mean(score, na.rm = TRUE),
          n_cells   = dplyr::n(),
          n_pos     = sum(score > 0),
          .groups   = "drop"
        )
      
      avg_scores_reactive(avg_scores)
      
      # -----------------------------
      # Cell type selector
      # -----------------------------
      # -----------------------------
      # Cell type selector + download button
      # -----------------------------
      output$celltype_selector_ui <- renderUI({
        req(avg_scores_reactive())
        avg_scores <- avg_scores_reactive()
        
        tagList(
          selectInput(
            session$ns("selected_annotation"),
            "Select cell type to visualize:",
            choices = sort(unique(avg_scores$annotation))
          ),
          br(),
          downloadButton(
            session$ns("download_table"),
            "Download Table"
          )
        )
      })
      
      
      # -----------------------------
      # Table output
      # -----------------------------
      output$metadata_table <- renderTable({
        avg_scores
      })
      
      # -----------------------------
      # Download button
      # -----------------------------
      output$download_table <- downloadHandler(
        filename = function() { paste0("UCell_scores_", Sys.Date(), ".csv") },
        content = function(file) {
          write.csv(avg_scores, file, row.names = FALSE)
        }
      )
      
    })
    
    # -----------------------------
    # Single bar plot for selected cell type
    # -----------------------------
    output$metadata_barplot <- renderPlot({
      req(avg_scores_reactive())
      
      df <- avg_scores_reactive()
      
      # Filter by selected cell type if any
      if (!is.null(input$selected_annotation)) {
        df <- df %>% dplyr::filter(annotation == input$selected_annotation)
      }
      
      validate(need(nrow(df) > 0, "No data to plot"))
      
      ggplot2::ggplot(df, ggplot2::aes(x = origin, y = avg_score, fill = origin)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "Dataset (origin.ident)",
          y = paste("Average UCell score")
        )
    })
    
  })
}
