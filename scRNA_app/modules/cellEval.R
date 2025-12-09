evalUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    fluidRow(
      column(
        4,
        selectInput(
          ns("species"),
          "Tissue origin:",
          choices = c("Human", "Mouse")
        )
      ),
      column(
        4,
        selectInput(
          ns("modality"),
          "Data type:",
          choices = c("Single-cell", "Single-nucleus")
        )
      )
    ),
    
    hr(),
    uiOutput(ns("file_status")),
    
    fluidRow(
      column(
        4,
        selectInput(ns("selected_celltype"), "Choose cell type:", choices = NULL),
        actionButton(ns("refresh_celltypes"), "Refresh cell types")
      ),
      column(
        8,
        downloadButton(ns("download_summary"), "Download gene % table")
      )
    ),
    
    br(),
    
    fluidRow(
      column(6, plotOutput(ns("gene_bar_plot"), height = "420px")),
      column(6, DT::DTOutput(ns("marker_count_table")))
    ),
    
    br(),
    h4("All gene % (long table)"),
    DT::DTOutput(ns("cell_table"))
  )
}


evalServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # ---- Marker file path ----
    marker_file <- reactive({
      sp <- input$species
      md <- input$modality
      req(sp, md)
      
      sp <- tolower(sp)
      md <- ifelse(md == "Single-cell", "sc", "sn")
      
      file.path(".", paste0(sp, " ", md, ".xlsx"))
    })
    
    # ---- Load marker table ----
    marker_df <- reactive({
      file <- marker_file()
      if (!file.exists(file)) return(NULL)
      
      ext <- tolower(tools::file_ext(file))
      
      df <- switch(
        ext,
        csv = try(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE),
        xlsx = try(readxl::read_excel(file), silent = TRUE),
        stop("Unsupported file type: ", ext)
      )
      
      if (inherits(df, "try-error") || ncol(df) < 2) return(NULL)
      
      df <- as.data.frame(df)
      df[] <- lapply(df, function(x) trimws(as.character(x)))
      
      colnames(df)[1:2] <- c("Tissue", "CellType")
      if (ncol(df) > 2)
        colnames(df)[3:ncol(df)] <- paste0("Marker", 1:(ncol(df) - 2))
      
      df
    })
    
    # ---- Percent expressing long table ----
    long_pct_table <- reactive({
      obj <- app_data$dataset %||% app_data$seurat_obj
      df  <- marker_df()
      if (is.null(df) || is.null(obj)) return(NULL)
      
      if (!"cell_ann" %in% colnames(obj@meta.data)) return(NULL)
      
      annotations <- obj@meta.data$cell_ann
      
      # Keep only cell types present in the dataset
      valid_rows <- df$CellType %in% unique(annotations)
      df <- df[valid_rows, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      
      out <- list()
      
      for (i in seq_len(nrow(df))) {
        
        cell_type <- df$CellType[i]
        markers <- as.character(unlist(df[i, 3:ncol(df)]))
        markers <- markers[nzchar(markers) & !is.na(markers)]
        if (length(markers) == 0) next
        
        present <- markers[markers %in% rownames(obj)]
        if (length(present) == 0) next
        
        cells_of_type <- rownames(obj@meta.data)[annotations == cell_type]
        if (length(cells_of_type) == 0) next
        
        expr_df <- tryCatch(
          FetchData(obj, vars = present, cells = cells_of_type),
          error = function(e) NULL
        )
        if (is.null(expr_df) || ncol(expr_df) == 0) next
        
        for (g in present) {
          pct <- mean(expr_df[[g]] > 0) * 100
          out[[length(out) + 1]] <- data.frame(
            CellType = cell_type,
            Gene = g,
            PercentExpressing = round(pct, 1),
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(out) == 0) return(NULL)
      do.call(rbind, out)
    })
    
    # ---- Available cell types ----
    matched_celltypes <- reactive({
      obj <- app_data$dataset %||% app_data$seurat_obj
      tbl <- marker_df()
      
      if (is.null(tbl) || is.null(obj)) return(character(0))
      if (!"cell_ann" %in% colnames(obj@meta.data)) return(character(0))
      
      intersect(unique(tbl$CellType), unique(obj@meta.data$cell_ann))
    })
    
    observe({
      types <- matched_celltypes()
      updateSelectInput(session, "selected_celltype",
                        choices = types,
                        selected = types[1])
    })
    
    observeEvent(input$refresh_celltypes, {
      types <- matched_celltypes()
      updateSelectInput(session, "selected_celltype",
                        choices = types,
                        selected = types[1])
    })
    
    # ---- Bar plot ----
    gene_bar_df <- reactive({
      req(input$selected_celltype)
      
      tbl <- long_pct_table()
      if (is.null(tbl)) return(NULL)
      
      df <- tbl[tbl$CellType == input$selected_celltype, , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      
      df$Gene <- factor(df$Gene, levels = df$Gene[order(df$PercentExpressing, decreasing = TRUE)])
      df
    })
    
    output$gene_bar_plot <- renderPlot({
      df <- gene_bar_df()
      
      if (is.null(df)) {
        ggplot() + theme_void() +
          ggtitle("No genes / data for selected cell type")
      } else {
        ggplot(df, aes(x = Gene, y = PercentExpressing)) +
          geom_col(fill = "#2C7BB6") +
          coord_flip() +
          ylab("% cells expressing") +
          xlab("") +
          ggtitle(sprintf("Percent expressing - %s", input$selected_celltype)) +
          theme_minimal()
      }
    })
    
    # ---- Marker-count distribution ----
    marker_count_distribution <- reactive({
      obj <- app_data$dataset %||% app_data$seurat_obj
      req(input$selected_celltype)
      
      df_markers <- marker_df()
      if (is.null(df_markers) || is.null(obj)) return(NULL)
      if (!"cell_ann" %in% colnames(obj@meta.data)) return(NULL)
      
      row <- df_markers[df_markers$CellType == input$selected_celltype, , drop = FALSE]
      
      markers <- as.character(unlist(row[1, 3:ncol(row)]))
      markers <- markers[nzchar(markers) & !is.na(markers)]
      if (length(markers) == 0) return(NULL)
      
      markers_present <- markers[markers %in% rownames(obj)]
      
      cells_of_type <- rownames(obj@meta.data)[
        obj@meta.data$cell_ann == input$selected_celltype
      ]
      
      expr_df <- tryCatch(
        FetchData(obj, vars = markers_present, cells = cells_of_type),
        error = function(e) NULL
      )
      if (is.null(expr_df)) return(NULL)
      
      counts <- rowSums(expr_df > 0)
      total_cells <- length(counts)
      N <- length(markers)
      
      freq_table <- as.data.frame(table(factor(counts, levels = 0:N)))
      colnames(freq_table) <- c("NumMarkersExpressed", "CellCount")
      
      freq_table$NumMarkersExpressed <- as.integer(as.character(freq_table$NumMarkersExpressed))
      freq_table$PercentCells <- round(100 * freq_table$CellCount / total_cells, 1)
      freq_table$FractionOfMarkers <- paste0(freq_table$NumMarkersExpressed, "/", N)
      
      freq_table
    })
    
    output$marker_count_table <- DT::renderDT({
      tbl <- marker_count_distribution()
      
      if (is.null(tbl)) {
        DT::datatable(
          data.frame(message = "No marker-count data for selected cell type."),
          options = list(dom = "t")
        )
      } else {
        DT::datatable(tbl, options = list(pageLength = 15), rownames = FALSE)
      }
    })
    
    # ---- Full long table ----
    output$cell_table <- DT::renderDT({
      tbl <- long_pct_table()
      
      if (is.null(tbl)) {
        DT::datatable(
          data.frame(message = "No matching cell types or genes."),
          options = list(dom = "t")
        )
      } else {
        DT::datatable(tbl, options = list(pageLength = 25), rownames = FALSE)
      }
    })
    
    # ---- Download ----
    output$download_summary <- downloadHandler(
      filename = function() paste0("marker_gene_percent_", Sys.Date(), ".csv"),
      content = function(file) {
        tbl <- long_pct_table()
        
        if (is.null(tbl)) {
          write.csv(data.frame(message = "No data"), file, row.names = FALSE)
        } else {
          write.csv(tbl, file, row.names = FALSE)
        }
      }
    )
    
    # ---- Status panel ----
    output$file_status <- renderUI({
      file <- marker_file()
      obj  <- app_data$dataset %||% app_data$seurat_obj
      
      if (!file.exists(file)) {
        return(
          wellPanel(
            h4("❌ Marker file not found"),
            p(sprintf("Expected: %s", file))
          )
        )
      }
      
      if (is.null(marker_df())) {
        return(
          wellPanel(
            h4("⚠️ File format problem"),
            p("Ensure file has a CellType column + marker genes.")
          )
        )
      }
      
      if (is.null(obj)) {
        return(
          wellPanel(
            h4("⚠️ Seurat object missing"),
            p("Load data first in the Upload step.")
          )
        )
      }
      
      if (!"cell_ann" %in% colnames(obj@meta.data)) {
        return(
          wellPanel(
            h4("⚠️ cell_ann missing"),
            p("Run annotation to create 'cell_ann' metadata.")
          )
        )
      }
      
      wellPanel(
        h4("✅ Marker file loaded"),
        p(sprintf("File: %s", file)),
        p(sprintf("Matching cell types: %d", length(matched_celltypes())))
      )
    })
  })
}
