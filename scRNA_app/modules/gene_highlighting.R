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
          actionButton(ns("add_filter"),    "➕ Add filter"),
          actionButton(ns("clear_filters"), "🧹 Clear filters", class = "btn-link")
        ),
        br(),
        uiOutput(ns("filters_ui")),   # dynamic filter blocks
        
        tags$hr(),
        uiOutput(ns("color_controls_ui")),  # color pickers (single / co / per-gene)
        
        actionButton(ns("plot_btn"), "Plot Feature"),
        br(), br(),
        downloadButton(ns("download_plot"),     "Download Plot"),
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
    
    # ---- small helpers ----
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    .active_obj <- reactive({
      if (!is.null(app_data$dataset)) app_data$dataset else app_data$seurat_obj
    })
    selected_genes <- reactive({
      gs <- unlist(strsplit(input$gene_name %||% "", "\\s+"))
      gs[nzchar(gs)]
    })
    
    # state
    plot_obj              <- reactiveVal(NULL)
    gene_used             <- reactiveVal(NULL)
    summary_data          <- reactiveVal(NULL)
    summary_by_origin     <- reactiveVal(NULL)
    summary_by_origin_ann <- reactiveVal(NULL)
    
    # ---------- Dynamic multi-filter system ----------
    filters_ids <- reactiveVal(character(0))
    next_id     <- reactiveVal(1)
    
    observeEvent(.active_obj(), {
      if (length(filters_ids()) == 0) {
        filters_ids(as.character(next_id()))
        next_id(next_id() + 1)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$add_filter, {
      id_new <- as.character(next_id())
      filters_ids(c(filters_ids(), id_new))
      next_id(next_id() + 1)
    })
    
    observeEvent(input$clear_filters, {
      filters_ids(character(0))
    })
    
    # Build the UI for all current filter blocks
    output$filters_ui <- renderUI({
      obj <- .active_obj()
      if (is.null(obj)) return(NULL)
      meta_cols <- colnames(obj@meta.data)
      
      tagList(lapply(filters_ids(), function(fid) {
        nsf <- function(x) session$ns(paste0(x, "_", fid))
        wellPanel(
          fluidRow(
            column(5,
                   selectInput(
                     nsf("field"),
                     "Metadata column:",
                     choices = meta_cols,
                     selected = meta_cols[1] %||% NULL
                   )
            ),
            column(7,
                   uiOutput(nsf("values_ui"))
            )
          )
        )
      }))
    })
    
    # Populate values per filter
    observe({
      obj <- .active_obj(); req(obj)
      meta_cols <- colnames(obj@meta.data)
      ids <- filters_ids()
      if (!length(ids)) return()
      
      lapply(ids, function(fid) {
        nsf <- function(x) session$ns(paste0(x, "_", fid))
        field_input <- reactive(input[[paste0("field_", fid)]])
        observeEvent(field_input(), {
          fld <- field_input()
          if (is.null(fld) || !(fld %in% meta_cols)) {
            output[[paste0("values_ui_", fid)]] <- renderUI(NULL)
            return()
          }
          vals <- sort(unique(as.character(obj@meta.data[[fld]])))
          output[[paste0("values_ui_", fid)]] <- renderUI({
            checkboxGroupInput(
              nsf("values"),
              label = "Values:",
              choices = vals,
              selected = vals  # select all by default
            )
          })
        }, ignoreInit = FALSE)
      })
    })
    
    # Cells that pass all filters
    filtered_cells <- reactive({
      obj <- .active_obj(); req(obj)
      if (!length(filters_ids())) return(colnames(obj))
      cells <- colnames(obj)
      for (fid in filters_ids()) {
        fld  <- input[[paste0("field_", fid)]]
        vals <- input[[paste0("values_", fid)]]
        if (is.null(fld) || is.null(vals) || !length(vals)) next
        hits <- which(as.character(obj@meta.data[[fld]]) %in% vals)
        cells <- intersect(cells, rownames(obj@meta.data)[hits])
      }
      cells
    })
    
    # ---------- Color controls ----------
    output$color_controls_ui <- renderUI({
      req(input$gene_name)
      genes <- selected_genes()
      if (!length(genes)) return(NULL)
      
      if (identical(input$co, "TRUE")) {
        colourpicker::colourInput(session$ns("co_color"), "Pick color for co-expressing cells:", value = "#0000FF")
      } else if (length(genes) == 1) {
        colourpicker::colourInput(session$ns("single_color"), "Pick color for expressing cells:", value = "#FF0000")
      } else {
        tagList(
          lapply(genes, function(g) {
            colourpicker::colourInput(session$ns(paste0("color_", g)), paste("Color for", g), value = "#0000FF")
          }),
          helpText("Non-expressing cells are always shown in gray.")
        )
      }
    })
    
    # ---------- Plot + Tables ----------
    observeEvent(input$plot_btn, {
      obj_full <- .active_obj(); req(obj_full)
      genes <- selected_genes()
      gene_used(genes)
      
      if (!length(genes)) { showNotification("❌ Please provide at least one gene.", type = "error"); return() }
      missing_genes <- setdiff(genes, rownames(obj_full))
      if (length(missing_genes) > 0) {
        showNotification(paste("❌ Not found:", paste(missing_genes, collapse = ", ")), type = "error")
        plot_obj(NULL); summary_data(NULL); summary_by_origin(NULL); summary_by_origin_ann(NULL)
        return()
      }
      
      cells_keep <- filtered_cells()
      if (!length(cells_keep)) {
        showNotification("⚠️ No cells match the selected filters.", type = "warning")
        plot_obj(NULL); summary_data(NULL); summary_by_origin(NULL); summary_by_origin_ann(NULL)
        return()
      }
      
      # Ensure UMAP exists
      if (is.null(obj_full@reductions$umap)) {
        showNotification("ℹ️ UMAP not found — running UMAP on full object…", type = "message")
        obj_full <- RunUMAP(obj_full, dims = 1:(app_data$num_pcs %||% 10))
        if (!is.null(app_data$dataset)) app_data$dataset <- obj_full else app_data$seurat_obj <- obj_full
      }
      
      umap_df <- data.frame(
        cell = cells_keep,
        Embeddings(obj_full, "umap")[cells_keep, , drop = FALSE]
      )
      colnames(umap_df) <- c("cell", "UMAP_1", "UMAP_2")
      
      if (length(genes) == 1 && input$co == "FALSE") {
        # single gene
        expr_vec <- FetchData(obj_full, vars = genes[1], cells = cells_keep)[, 1]
        df <- cbind(umap_df, expr = expr_vec)
        chosen_color <- input$single_color %||% "#FF0000"
        
        p <- ggplot() +
          geom_point(data = df, aes(UMAP_1, UMAP_2), color = "lightgray", size = 0.6) +
          geom_point(
            data = df[df$expr > 0, , drop = FALSE],
            aes(UMAP_1, UMAP_2, alpha = expr),
            color = chosen_color, size = 0.6
          ) +
          scale_alpha_continuous(range = c(0.3, 1)) +
          ggtitle(paste("Expression of:", genes[1])) +
          theme_minimal() + theme(legend.position = "none")
        plot_obj(p)
        
      } else if (identical(input$co, "TRUE")) {
        # co-expression (all genes > 0)
        expr_mat <- FetchData(obj_full, vars = genes, cells = cells_keep)
        co_cells <- rownames(expr_mat)[apply(expr_mat > 0, 1, all)]
        umap_df$highlight <- ifelse(umap_df$cell %in% co_cells, "Yes", "No")
        chosen_color <- input$co_color %||% "#0000FF"
        
        p <- ggplot() +
          geom_point(data = umap_df[umap_df$highlight == "No", ], aes(UMAP_1, UMAP_2),
                     color = "lightgray", size = 0.6) +
          geom_point(data = umap_df[umap_df$highlight == "Yes", ], aes(UMAP_1, UMAP_2),
                     color = chosen_color, size = 0.6) +
          ggtitle(paste("Cells co-expressing:", paste(genes, collapse = ", "))) +
          theme_minimal()
        plot_obj(p)
        
      } else {
        # multi-gene overlay
        expr_df <- FetchData(obj_full, vars = genes, cells = cells_keep)
        expr_df$cell <- rownames(expr_df)
        combined <- merge(umap_df, expr_df, by = "cell")
        long_df <- reshape2::melt(
          combined,
          id.vars = c("cell", "UMAP_1", "UMAP_2"),
          variable.name = "gene", value.name = "expression"
        )
        expressed_df <- long_df[long_df$expression > 0, ]
        
        gene_colors <- sapply(genes, function(g) input[[paste0("color_", g)]] %||% "#0000FF")
        color_map <- setNames(gene_colors, genes)
        
        p <- ggplot() +
          geom_point(data = combined, aes(UMAP_1, UMAP_2), color = "lightgray", size = 0.6) +
          geom_point(
            data = expressed_df,
            aes(UMAP_1, UMAP_2, color = gene, alpha = expression),
            size = 0.6
          ) +
          scale_color_manual(values = color_map) +
          scale_alpha_continuous(range = c(0.3, 1)) +
          ggtitle(paste("Expression of:", paste(genes, collapse = ", "))) +
          theme_minimal() + theme(legend.title = element_blank())
        plot_obj(p)
      }
      
      # ---- Summary tables (on filtered cells) ----
      expr_mat <- FetchData(obj_full, vars = genes, cells = cells_keep)
      meta     <- obj_full@meta.data[cells_keep, , drop = FALSE]
      meta$cell <- rownames(meta)
      expr_mat$cell <- rownames(expr_mat)
      combined <- merge(meta, expr_mat, by = "cell")
      
      if (!"cell_ann" %in% colnames(combined)) {
        combined$cell_ann <- as.character(combined$seurat_clusters)
      }
      
      origin_col <- if ("orig.ident" %in% colnames(combined)) {
        "orig.ident"
      } else if ("sample" %in% colnames(combined)) {
        "sample"
      } else {
        "seurat_clusters"
      }
      combined$.origin <- as.character(combined[[origin_col]])
      combined$any_expr <- apply(combined[, genes, drop = FALSE] > 0, 1, any)
      combined$co_expr  <- apply(combined[, genes, drop = FALSE] > 0, 1, all)
      
      summary_data(
        combined %>%
          dplyr::group_by(cell_ann) %>%
          dplyr::summarise(
            total_cells      = dplyr::n(),
            any_expr_cells   = sum(any_expr),
            co_expr_cells    = sum(co_expr),
            any_expr_percent = round(100 * mean(any_expr), 1),
            co_expr_percent  = round(100 * mean(co_expr), 1),
            dplyr::across(dplyr::all_of(genes), ~mean(.x[.x > 0]), .names = "avg_{.col}_expressing")
          ) %>%
          dplyr::arrange(dplyr::desc(any_expr_percent))
      )
      
      summary_by_origin(
        combined %>%
          dplyr::group_by(.origin) %>%
          dplyr::summarise(
            total_cells      = dplyr::n(),
            any_expr_cells   = sum(any_expr),
            co_expr_cells    = sum(co_expr),
            any_expr_percent = round(100 * mean(any_expr), 1),
            co_expr_percent  = round(100 * mean(co_expr), 1),
            dplyr::across(dplyr::all_of(genes), ~mean(.x[.x > 0]), .names = "avg_{.col}_expressing")
          ) %>%
          dplyr::arrange(dplyr::desc(any_expr_percent))
      )
      
      summary_by_origin_ann(
        combined %>%
          dplyr::group_by(.origin, cell_ann) %>%
          dplyr::summarise(
            total_cells      = dplyr::n(),
            any_expr_cells   = sum(any_expr),
            co_expr_cells    = sum(co_expr),
            any_expr_percent = round(100 * mean(any_expr), 1),
            co_expr_percent  = round(100 * mean(co_expr), 1),
            dplyr::across(dplyr::all_of(genes), ~mean(.x[.x > 0]), .names = "avg_{.col}_expressing"),
            .groups = "drop"
          ) %>%
          dplyr::arrange(.origin, dplyr::desc(any_expr_percent))
      )
    })
    
    # ---------- Helper for raw expression slice (for download) ----------
    expr_slice <- reactive({
      obj <- .active_obj(); req(obj)
      gs <- selected_genes(); req(length(gs) > 0)
      keep <- filtered_cells(); req(length(keep) > 0)
      # FetchData returns a data.frame with columns = genes and rows = cells
      FetchData(obj, vars = gs, cells = keep)
    })
    
    # ---------- Outputs ----------
    output$feature_plot <- renderPlot({ req(plot_obj()); plot_obj() })
    output$expression_summary <- DT::renderDT({
      req(summary_data())
      DT::datatable(summary_data(), options = list(pageLength = 10), rownames = FALSE)
    })
    output$expression_summary_by_origin <- DT::renderDT({
      req(summary_by_origin())
      DT::datatable(summary_by_origin(), options = list(pageLength = 10), rownames = FALSE)
    })
    output$expression_summary_by_origin_ann <- DT::renderDT({
      req(summary_by_origin_ann())
      DT::datatable(summary_by_origin_ann(), options = list(pageLength = 12), rownames = FALSE)
    })
    
    # Summary box (counts + averages among expressing cells)
    output$cell_summary_box <- renderUI({
      obj <- .active_obj(); req(obj, gene_used())
      genes <- gene_used()
      
      cells_keep <- filtered_cells()
      if (!length(cells_keep)) {
        return(wellPanel(HTML("<i>No cells match the current filters.</i>")))
      }
      
      expr_mat <- FetchData(obj, vars = genes, cells = cells_keep)
      total_cells <- nrow(expr_mat)
      any_expr    <- apply(expr_mat > 0, 1, any)
      expressing_cells <- sum(any_expr)
      percent_expressing <- if (total_cells > 0) round(100 * expressing_cells / total_cells, 2) else 0
      
      per_gene_avgs <- vapply(genes, function(g) {
        v <- expr_mat[[g]]
        if (is.null(v)) return(NA_real_)
        pos <- v[v > 0]
        if (length(pos) == 0) return(NA_real_)
        round(mean(pos), 3)
      }, numeric(1))
      
      avg_list_html <- paste0(
        "<ul style='margin:6px 0 0 18px;padding:0;'>",
        paste(sprintf("<li><b>%s</b>: %s", htmltools::htmlEscape(names(per_gene_avgs)),
                      ifelse(is.na(per_gene_avgs), "—", format(per_gene_avgs, trim = TRUE))),
              collapse = "</li>"),
        "</li></ul>"
      )
      
      wellPanel(
        HTML(sprintf(
          "<b>Total (filtered) cells:</b> %d<br>
           <b>Cells expressing any selected gene:</b> %d (%.2f%% of filtered cells)<br>
           <b>Average expression among expressing cells (per gene):</b>%s",
          total_cells, expressing_cells, percent_expressing, avg_list_html
        ))
      )
    })
    
    # ---------- Downloads ----------
    # Plot
    output$download_plot <- downloadHandler(
      filename = function() paste0("multi_gene_expression_", paste(gene_used(), collapse = "_"), "_", Sys.Date(), ".png"),
      content = function(file) {
        req(plot_obj())
        ggsave(file, plot = plot_obj(), width = 10, height = 8, dpi = 300)
      }
    )
    
    # Raw expression: filename includes genes + active filters
    output$download_raw_expr <- downloadHandler(
      filename = function() {
        gs <- paste(selected_genes(), collapse = "_"); if (!nzchar(gs)) gs <- "genes"
        active_filters <- c()
        for (fid in filters_ids()) {
          fld  <- input[[paste0("field_", fid)]]
          vals <- input[[paste0("values_", fid)]]
          if (!is.null(fld) && length(vals) > 0) {
            filter_tag <- paste0(fld, "=", paste(vals, collapse = "-"))
            filter_tag <- gsub("[^A-Za-z0-9_-]", "_", filter_tag)  # sanitize
            active_filters <- c(active_filters, filter_tag)
          }
        }
        filter_str <- if (length(active_filters)) paste(active_filters, collapse = "_") else "all_cells"
        paste0("expr_", gs, "_", filter_str, "_", Sys.Date(), ".csv.gz")
      },
      content = function(file) {
        df <- expr_slice()  # cells x genes
        df_out <- cbind(cell = rownames(df), df)
        con <- gzfile(file, open = "wb")
        on.exit(close(con), add = TRUE)
        utils::write.csv(df_out, con, row.names = FALSE)
      }
    )
  })
}