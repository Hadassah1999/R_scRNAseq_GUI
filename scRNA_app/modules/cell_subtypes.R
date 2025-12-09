subtypeUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      radioButtons(
        ns("species"),
        "Select species:",
        choices = c("Human", "Mouse"),
        selected = "Human"
      ),
      actionButton(ns("find_subtypes"), "Find Subtypes"),
      uiOutput(ns("select_celltype_ui"))
    ),
    mainPanel(
      plotOutput(ns("subtype_barplot")),
      tableOutput(ns("subtype_table"))
    )
  )
}

subtypeServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # -----------------------------
    # LOAD EXCEL
    # -----------------------------
    subtype_df <- eventReactive(input$find_subtypes, {
      req(input$species)
      file_path <- file.path("Cell_Subtypes", paste0(input$species, "_cells.xlsx"))
      print(paste("Looking for file:", file_path))
      req(file.exists(file_path))
      df <- readxl::read_excel(file_path)
      print(paste("Excel file loaded. Rows:", nrow(df), "Cols:", ncol(df)))
      df
    })
    
    # -----------------------------
    # ANNOTATION
    # -----------------------------
    subtype_results <- eventReactive(input$find_subtypes, {
      notif_id <- showNotification("Calculating subtypes...", duration = NULL, type = "message")
      on.exit(removeNotification(notif_id))
      
      obj <- app_data$dataset %||% app_data$seurat_obj
      req(is(obj, "Seurat"))
      req(subtype_df())
      
      meta <- obj@meta.data
      df <- subtype_df()
      
      # Keep only cell types present in data
      df <- df[df$`Cell Type` %in% meta$cell_ann, ]
      
      # Build marker list
      subtype_list <- df %>%
        tidyr::pivot_longer(
          cols = starts_with("Marker Gene"),
          names_to = "marker_number",
          values_to = "gene"
        ) %>%
        filter(!is.na(gene)) %>%
        group_by(`Cell Type`, Subtype) %>%
        summarize(markers = list(unique(gene)), .groups = "drop")
      
      print("Marker list prepared:")
      print(subtype_list)
      
      # Prepare results dataframe
      results <- data.frame(
        cell = rownames(meta),
        cell_ann = meta$cell_ann,
        orig.ident = meta$orig.ident,
        subtype = NA,
        stringsAsFactors = FALSE
      )
      
      # Loop over each cell type
      for (cell_type in unique(meta$cell_ann)) {
        print(paste("Annotating cell type:", cell_type))
        
        cells_of_type <- rownames(meta)[meta$cell_ann == cell_type]
        if (length(cells_of_type) == 0) next
        
        subtypes_for_type <- subtype_list %>% filter(`Cell Type` == cell_type)
        if (nrow(subtypes_for_type) == 0) next
        
        score_list <- list()
        
        for (i in seq_len(nrow(subtypes_for_type))) {
          marker_genes <- subtypes_for_type$markers[[i]]
          genes_present <- intersect(marker_genes, rownames(obj))
          if (length(genes_present) == 0) next
          
          expr_mat <- tryCatch(FetchData(obj, vars = genes_present, cells = cells_of_type), error = function(e) NULL)
          if (is.null(expr_mat) || nrow(expr_mat) == 0) next
          
          score_list[[subtypes_for_type$Subtype[i]]] <- rowSums(expr_mat > 0) / length(marker_genes)
          
        }
        
        if (length(score_list) == 0) next
        
        score_mat <- do.call(cbind, score_list)
        
        # ✅ Assign max-score subtype, but NA if no markers detected
        assigned <- apply(score_mat, 1, function(x) {
          if (all(x == 0)) return(NA)
          names(x)[which.max(x)]
        })
        
        results$subtype[match(names(assigned), results$cell)] <- assigned
      }
      
      # Assign Unclassified for cells with NA
      results$subtype[is.na(results$subtype)] <- "Unclassified"
      
      print("Annotation finished.")
      head(results)
      
      results
    })
    
    # -----------------------------
    # UI: VALID CELL TYPES
    # -----------------------------
    matched_celltypes <- reactive({
      req(subtype_results())
      valid <- subtype_results() %>%
        dplyr::filter(!is.na(subtype)) %>%
        dplyr::pull(cell_ann) %>%
        unique()
      
      if (length(valid) == 0) return(NULL)
      valid
    })
    
    output$select_celltype_ui <- renderUI({
      req(matched_celltypes())
      tagList(
        br(),
        selectInput(ns("selected_celltype"), "Select cell type:", choices = matched_celltypes()),
        p("Notice: If all of the cells are labeled 'Unclassified', you may have selected the wrong species.")
      )
    })
    
    # -----------------------------
    # BARPLOT
    # -----------------------------
    output$subtype_barplot <- renderPlot({
      req(subtype_results())
      req(input$selected_celltype)
      
      df <- subtype_results() %>%
        dplyr::filter(cell_ann == input$selected_celltype)
      
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No subtypes detected for this cell type")
        return(NULL)
      }
      
      summary_df <- df %>%
        dplyr::group_by(subtype) %>%
        dplyr::summarize(n_cells = n(), .groups = "drop")
      
      ggplot2::ggplot(summary_df, ggplot2::aes(x = subtype, y = n_cells, fill = subtype)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = paste("Subtypes of", input$selected_celltype), y = "Cells") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })
    
    # -----------------------------
    # TABLE
    # -----------------------------
    output$subtype_table <- renderTable({
      req(subtype_results())
      req(input$selected_celltype)
      
      df <- subtype_results() %>%
        dplyr::filter(cell_ann == input$selected_celltype)   # include all subtypes
      
      if (nrow(df) == 0) return(NULL)
      
      summary_long <- df %>%
        dplyr::group_by(subtype, orig.ident) %>%
        dplyr::summarize(n = n(), .groups = "drop")
      
      summary_wide <- summary_long %>%
        tidyr::pivot_wider(
          names_from = orig.ident,
          values_from = n,
          values_fill = 0
        )
      
      # Add TOTAL row
      total_row <- summary_wide %>%
        dplyr::summarize(
          subtype = "**TOTAL**",
          dplyr::across(where(is.numeric), sum)
        )
      
      final_table <- dplyr::bind_rows(summary_wide, total_row)
      final_table
    })
    
  })
}
