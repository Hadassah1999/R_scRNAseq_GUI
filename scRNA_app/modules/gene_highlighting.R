geneHighlightingUI <- function(id) {
  ns <- NS(id)
  tagList(
    titlePanel("Gene Expression Viewer"),
    sidebarLayout(
      sidebarPanel(
        textInput(ns("gene_name"), "Enter gene names (space-separated, case-sensitive):"),
        radioButtons(
          ns("co"),
          "Show only co-expressed genes?",
          choices = c("yes" = "TRUE", "no" = "FALSE"),
          selected = "FALSE"
        ),
        tags$hr(),
        h4("Filters"),
        tagList(
          actionButton(ns("add_filter"), "➕ Add filter"),
          actionButton(ns("clear_filters"), "🧹 Clear filters", class = "btn-link")
        ),
        br(),
        uiOutput(ns("filters_ui")),
        tags$hr(),
        uiOutput(ns("color_controls_ui")),
        checkboxInput(ns("show_clusters"), "Show cluster numbers", value = FALSE),
        checkboxInput(ns("show_annotations"), "Show annotation labels", value = FALSE),
        actionButton(ns("plot_btn"), "Plot Feature"),
        br(), br(),
        downloadButton(ns("download_plot"), "Download Plot"),
        br(), br(),
        downloadButton(ns("download_raw_expr"), "Download Raw Expression (CSV.GZ)")
      ),
      mainPanel(
        h4("FeaturePlot"),
        withSpinner(plotOutput(ns("feature_plot")), type = 6),
        br(),
        uiOutput(ns("cell_summary_box")),
        br(),
        h4("Summary Table (by annotation)"),
        DT::DTOutput(ns("expression_summary")),
        br(),
        h4("Summary Table (by origin)"),
        DT::DTOutput(ns("expression_summary_by_origin")),
        br(),
        h4("Summary Table (by origin × annotation)"),
        DT::DTOutput(ns("expression_summary_by_origin_ann"))
      )
    )
  )
}

geneHighlightingServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    
    # ---- Reactives ----
    .active_obj <- reactive({ app_data$dataset %||% app_data$seurat_obj })
    
    selected_genes <- reactive({
      genes <- strsplit(input$gene_name %||% "", "\\s+")[[1]]
      genes <- trimws(genes)
      genes[genes != ""]
    })
    
    plot_obj <- reactiveVal(NULL)
    gene_used <- reactiveVal(NULL)
    summary_data <- reactiveVal(NULL)
    summary_by_origin <- reactiveVal(NULL)
    summary_by_origin_ann <- reactiveVal(NULL)
    
    # ---------- Dynamic filters ----------
    filters_ids <- reactiveVal(character(0))
    next_id <- reactiveVal(1)
    
    init_default_filter <- function() {
      filters_ids("1")
      next_id(2)
    }
    
    observeEvent(.active_obj(), {
      req(.active_obj())
      # clear previous filters and start fresh
      filters_ids(character(0))
      next_id(1)
      init_default_filter()
    }, ignoreInit = FALSE)
    

    
    observeEvent(input$add_filter, {
      id_new <- as.character(next_id())
      filters_ids(c(filters_ids(), id_new))
      next_id(next_id() + 1)
      
      # Set a default metadata column for the new filter
      obj <- .active_obj(); req(obj)
      meta_cols <- colnames(obj@meta.data)
      if (length(meta_cols)) {
        updateSelectInput(session, paste0("field_", id_new), selected = meta_cols[1])
      }
    })
    
    observeEvent(input$clear_filters, {
      init_default_filter()
    })
    
    # Build UI for dynamic filters
    output$filters_ui <- renderUI({
      obj <- .active_obj(); req(obj)
      meta_cols <- colnames(obj@meta.data)
      tagList(
        lapply(filters_ids(), function(fid) {
          nsf <- function(x) session$ns(paste0(x, "_", fid))
          wellPanel(
            fluidRow(
              column(
                5,
                selectInput(nsf("field"), "Metadata column:", choices = meta_cols, selected = meta_cols[1])
              ),
              column(7, uiOutput(nsf("values_ui")))
            )
          )
        })
      )
    })
    
    # Populate values for each filter safely
    observe({
      obj <- .active_obj(); req(obj)
      meta_cols <- colnames(obj@meta.data)
      
      lapply(filters_ids(), function(fid) {
        nsf <- function(x) session$ns(paste0(x, "_", fid))
        
        field_val <- input[[paste0("field_", fid)]]
        if (is.null(field_val) || !(field_val %in% meta_cols)) return(NULL)
        
        vals <- sort(unique(as.character(obj@meta.data[[field_val]])))
        
        output[[paste0("values_ui_", fid)]] <- renderUI({
          tagList(
            checkboxGroupInput(nsf("values"), "Values:", choices = vals, selected = vals),
            checkboxInput(nsf("select_all"), "Select All", value = TRUE)
          )
        })
        
        observeEvent(input[[paste0("select_all_", fid)]], {
          if (isTRUE(input[[paste0("select_all_", fid)]])) {
            updateCheckboxGroupInput(session, paste0("values_", fid), selected = vals)
          } else {
            updateCheckboxGroupInput(session, paste0("values_", fid), selected = character(0))
          }
        }, ignoreNULL = FALSE)
      })
    })
    
    # Filtered cells based on metadata filters
    filtered_cells <- reactive({
      obj <- .active_obj(); req(obj)
      cells <- colnames(obj)
      for (fid in filters_ids()) {
        fld  <- input[[paste0("field_", fid)]]
        vals <- input[[paste0("values_", fid)]]
        if (!is.null(fld) && length(vals)) {
          hits  <- which(as.character(obj@meta.data[[fld]]) %in% vals)
          cells <- intersect(cells, rownames(obj@meta.data)[hits])
        }
      }
      cells
    })
    
    # ---------- Color persistence ----------
    color_store <- reactiveVal(list())
    
    # ---- Random color helpers ----
    .golden_h <- reactiveVal(runif(1))  # random start per session
    next_distinct_hex <- function() {
      h <- (.golden_h() + 0.61803398875) %% 1
      .golden_h(h)
      grDevices::hcl(h * 360, 70, 55)   # readable HCL color
    }
    
    
    # Watch for user color input changes and store them persistently
    observe({
      cols <- color_store()
      # capture all color-related inputs
      new_cols <- reactiveValuesToList(input)
      new_cols <- new_cols[grepl("^color_|^co_color$|^single_color$", names(new_cols))]
      # merge updates into color_store
      color_store(modifyList(cols, new_cols))
    })
    
    # Track previous gene list to detect newly added genes
    .prev_genes <- reactiveVal(character(0))
    
    # When genes or 'co' change, if multi & non-co, assign random colors to new genes
    observeEvent(list(selected_genes(), input$co), {
      gs <- selected_genes()
      old <- .prev_genes()
      
      # Always update the tracker
      .prev_genes(gs)
      
      # Only act in multi-gene, non-co mode
      if (!identical(input$co, "FALSE") || length(gs) <= 1) return()
      
      added <- setdiff(gs, old)
      if (!length(added)) return()
      
      cols <- color_store()
      for (g in added) {
        key <- paste0("color_", g)
        if (is.null(cols[[key]]) || !nzchar(cols[[key]])) {
          cols[[key]] <- next_distinct_hex()
        }
      }
      color_store(cols)
    }, ignoreInit = TRUE)
    
    
    # ---------- Color controls ----------
    output$color_controls_ui <- renderUI({
      genes <- selected_genes(); req(length(genes))
      cols  <- color_store()
      
      if (identical(input$co, "TRUE")) {
        colourpicker::colourInput(
          session$ns("co_color"),
          "Pick color for co-expressing cells:",
          value = cols$co_color %||% "#0000FF"
        )
      } else if (length(genes) == 1) {
        colourpicker::colourInput(
          session$ns("single_color"),
          "Pick color for expressing cells:",
          value = cols$single_color %||% "#0000FF"
        )
      } else {
        tagList(
          lapply(genes, function(g) {
            colourpicker::colourInput(
              session$ns(paste0("color_", g)),
              paste("Color for", g),
              value = cols[[paste0("color_", g)]] %||% "#0000FF"
              
              
            )
          }),
          helpText("Non-expressing cells are always shown in gray.")
        )
      }
    })
    
    # ---------- Plot & Summaries ----------
    create_summary <- function(df, group_vars) {
      df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
          total_cells = dplyr::n(),
          any_expr_cells = sum(any_expr),
          dplyr::across(dplyr::all_of(selected_genes()), ~ sum(.x > 0), .names = "expr_cells_{.col}"),
          co_expr_cells = sum(co_expr),
          any_expr_percent = round(100 * mean(any_expr), 1),
          co_expr_percent = round(100 * mean(co_expr), 1),
          dplyr::across(dplyr::all_of(selected_genes()), ~ mean(.x[.x > 0]), .names = "avg_{.col}_expressing"),
          .groups = "drop"
        )
    }
    
    observeEvent(input$plot_btn, {
      obj <- .active_obj(); req(obj)
      
      # --- Require cell annotations before proceeding ---
      if (!"cell_ann" %in% colnames(obj@meta.data) ||
          all(is.na(obj$cell_ann)) ||
          all(trimws(as.character(obj$cell_ann)) == "")) {
        showNotification(
          "❌ Cell annotation not found. Please run the cell annotation module before using this feature.",
          type = "error",
          duration = 10
        )
        return(invisible(NULL))  # stop execution before plotting
      }
      
      genes <- selected_genes(); gene_used(genes)
      
      # If genes are provided, check if they exist
      if (length(genes)) {
        missing_genes <- setdiff(genes, rownames(obj))
        if (length(missing_genes))
          return(showNotification(paste("❌ Not found:", paste(missing_genes, collapse = ", ")), type = "error"))
      }
      
      cells_keep <- filtered_cells()
      if (!length(cells_keep))
        return(showNotification("⚠️ No cells match the selected filters.", type = "warning"))
      
      # Ensure UMAP exists
      if (!"umap" %in% names(obj@reductions)) {
        obj <- RunUMAP(obj, dims = 1:(app_data$num_pcs %||% 10), verbose = FALSE)
        # persist
        if (!is.null(app_data$dataset)) app_data$dataset <- obj else app_data$seurat_obj <- obj
      }
      
      umap_df <- data.frame(cell = cells_keep, Embeddings(obj, "umap")[cells_keep, , drop = FALSE])
      colnames(umap_df) <- c("cell", "UMAP_1", "UMAP_2")
      
      # Cluster and annotation info
      umap_df$cluster    <- as.character(obj$seurat_clusters[cells_keep])
      umap_df$annotation <- obj$cell_ann[cells_keep] %||% umap_df$cluster
      umap_df$annotation <- ifelse(is.na(umap_df$annotation), umap_df$cluster, umap_df$annotation)
      
      # Compute centroids for labels
      cluster_centroids <- umap_df %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
      
      annotation_centroids <- umap_df %>%
        dplyr::group_by(annotation) %>%
        dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
      
      # ---- Plot ----
      plot_obj({
        # If no genes selected → just regular gray UMAP
        if (!length(genes)) {
          p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
            geom_point(color = "lightgray", size = 0.6) +
            theme_minimal() +
            ggtitle("UMAP (no genes selected)")
          
        } else if (length(genes) == 1 && input$co == "FALSE") {
          # Single gene expression
          expr_vec <- FetchData(obj, vars = genes[1], cells = cells_keep)[, 1]
          df <- cbind(umap_df, expr = expr_vec)
          p <- ggplot() +
            geom_point(data = umap_df, aes(UMAP_1, UMAP_2), color = "lightgray", size = 0.6) +
            geom_point(
              data  = df[df$expr > 0, ],
              aes(UMAP_1, UMAP_2, alpha = expr),
              color = input$single_color %||% "#0000FF",
              size  = 0.6
            ) +
            scale_alpha_continuous(range = c(0.3, 1)) +
            theme_minimal() +
            ggtitle(paste("Expression of:", genes[1]))
          
        } else if (identical(input$co, "TRUE")) {
          # Co-expression
          expr_mat  <- FetchData(obj, vars = genes, cells = cells_keep)
          co_cells  <- rownames(expr_mat)[apply(expr_mat > 0, 1, all)]
          umap_df$highlight <- ifelse(umap_df$cell %in% co_cells, "Yes", "No")
          p <- ggplot() +
            geom_point(data = umap_df, aes(UMAP_1, UMAP_2), color = "lightgray", size = 0.6) +
            geom_point(
              data  = umap_df[umap_df$highlight == "Yes", ],
              aes(UMAP_1, UMAP_2),
              color = input$co_color %||% "#0000FF",
              size  = 0.6
            ) +
            theme_minimal() +
            ggtitle(paste("Co-expression of:", paste(genes, collapse = ", ")))
          
        } else {
          # Multiple genes (different colors)
          expr_df  <- FetchData(obj, vars = genes, cells = cells_keep)
          combined <- cbind(umap_df, expr_df)
          long_df  <- tidyr::pivot_longer(
            combined, cols = genes, names_to = "gene", values_to = "expression"
          )
          p <- ggplot() +
            geom_point(data = umap_df, aes(UMAP_1, UMAP_2), color = "lightgray", size = 0.6) +
            geom_point(
              data  = long_df[long_df$expression > 0, ],
              aes(UMAP_1, UMAP_2, color = gene, alpha = expression),
              size  = 0.6
            ) +
            scale_color_manual(
              values = setNames(
                sapply(genes, function(g) {
                  color_store()[[paste0("color_", g)]] %||% next_distinct_hex()
                }),
                genes
              )
            ) +
            scale_alpha_continuous(range = c(0.3, 1)) +
            theme_minimal() +
            ggtitle(paste("Expression of:", paste(genes, collapse = ", ")))
        }
        
        # ---- ADD CLUSTER / ANNOTATION LABELS TO ALL CASES ----
        if (isTRUE(input$show_clusters))
          p <- p + geom_label(
            data = cluster_centroids,
            aes(UMAP_1, UMAP_2, label = cluster),
            color = "black", fontface = "bold", size = 4
          )
        
        if (isTRUE(input$show_annotations))
          p <- p + geom_label(
            data = annotation_centroids,
            aes(UMAP_1, UMAP_2, label = annotation),
            color = "darkgreen", fontface = "bold", size = 3
          )
        
        p
      })
      
      
      # ---- Summaries (only if genes selected) ----
      if (length(genes)) {
        expr_mat <- FetchData(obj, vars = genes, cells = cells_keep)
        meta     <- obj@meta.data[cells_keep, , drop = FALSE]
        combined <- cbind(meta, expr_mat)
        combined$any_expr <- apply(combined[, genes, drop = FALSE] > 0, 1, any)
        combined$co_expr  <- apply(combined[, genes, drop = FALSE] > 0, 1, all)
        combined$cell_ann <- combined$cell_ann %||% as.character(combined$seurat_clusters)
        combined$.origin  <- combined$orig.ident %||% combined$sample %||% as.character(combined$seurat_clusters)
        
        summary_data(create_summary(combined, "cell_ann"))
        summary_by_origin(create_summary(combined, ".origin"))
        summary_by_origin_ann(create_summary(combined, c(".origin", "cell_ann")))
      }
    })
    
    # ---------- Cell summary ----------
    output$cell_summary_box <- renderUI({
      obj <- .active_obj(); req(gene_used())
      genes <- gene_used()
      cells_keep <- filtered_cells()
      if (!length(cells_keep)) return(wellPanel(HTML("<i>No cells match the current filters.</i>")))
      expr_mat <- FetchData(obj, vars = genes, cells = cells_keep)
      any_expr <- apply(expr_mat > 0, 1, any)
      per_gene_avgs <- purrr::map_dbl(genes, ~ mean(expr_mat[[.x]][expr_mat[[.x]] > 0], na.rm = TRUE))
      avg_list_html <- paste0(
        "<ul style='margin:6px 0 0 18px;'>",
        paste(sprintf("<li><b>%s</b>: %s", genes, format(per_gene_avgs, trim = TRUE)), collapse = "</li>"),
        "</li></ul>"
      )
      wellPanel(HTML(sprintf(
        "<b>Total (filtered) cells:</b> %d<br><b>Cells expressing any selected gene:</b> %d (%.2f%%)<br><b>Average expression per gene:</b>%s",
        nrow(expr_mat), sum(any_expr), 100 * mean(any_expr), avg_list_html
      )))
    })
    
    # ---------- Outputs ----------
    output$feature_plot <- renderPlot({ req(plot_obj()); plot_obj() })
    output$expression_summary <- DT::renderDT({
      req(summary_data()); DT::datatable(summary_data(), options = list(pageLength = 10), rownames = FALSE)
    })
    output$expression_summary_by_origin <- DT::renderDT({
      req(summary_by_origin()); DT::datatable(summary_by_origin(), options = list(pageLength = 10), rownames = FALSE)
    })
    output$expression_summary_by_origin_ann <- DT::renderDT({
      req(summary_by_origin_ann()); DT::datatable(summary_by_origin_ann(), options = list(pageLength = 12), rownames = FALSE)
    })
    
    # ---------- Downloads ----------
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("multi_gene_expression_", paste(gene_used(), collapse = "_"), "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(plot_obj())
        # Ensure the saved plot has a white background
        ggsave(
          filename = file,
          plot = plot_obj() + theme(plot.background = element_rect(fill = "white", color = NA)),
          width = 10,
          height = 8,
          dpi = 300,
          bg = "white"  # explicitly set background
        )
      }
    )
    
    expr_slice <- reactive({
      obj <- .active_obj(); req(length(selected_genes()) > 0)
      FetchData(obj, vars = selected_genes(), cells = filtered_cells())
    })
    
    output$download_raw_expr <- downloadHandler(
      filename = function() {
        gs <- paste(selected_genes(), collapse = "_"); if (!nzchar(gs)) gs <- "genes"
        filter_tags <- purrr::map_chr(filters_ids(), function(fid) {
          fld  <- input[[paste0("field_", fid)]]
          vals <- input[[paste0("values_", fid)]]
          if (!is.null(fld) && length(vals)) paste0(fld, "=", paste(vals, collapse = "-")) else NULL
        }) |> na.omit()
        paste0(
          "expr_", gs, "_",
          ifelse(length(filter_tags), paste(filter_tags, collapse = "_"), "all_cells"),
          "_", Sys.Date(), ".csv.gz"
        )
      },
      content = function(file) {
        df <- expr_slice()
        df <- cbind(cell = rownames(df), df)
        write.csv(df, gzfile(file), row.names = FALSE)
      }
    )
    
    # ---------- Reset handler ----------
    observeEvent(app_data$reset_gene_highlighting, {
      # Reset reactive values
      plot_obj(NULL)
      gene_used(NULL)
      summary_data(NULL)
      summary_by_origin(NULL)
      summary_by_origin_ann(NULL)
      color_store(list())
      
      # Reinitialize dynamic filters
      filters_ids(character(0))
      next_id(1)
      init_default_filter()
      
      # Reset inputs
      updateTextInput(session, "gene_name", value = "")
      updateRadioButtons(session, "co", selected = "FALSE")
      updateCheckboxInput(session, "show_clusters", value = FALSE)
      updateCheckboxInput(session, "show_annotations", value = FALSE)

      
    })
  })
}
