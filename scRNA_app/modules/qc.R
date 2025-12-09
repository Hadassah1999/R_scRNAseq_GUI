qcModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Add Sample Paths"),
    div(
      style = "display: flex; flex-direction: column; gap: 6px; max-width: 500px;",
      textInput(ns("new_path"), "Enter path to .h5 or 10X folder (matrix dir):"),
      textInput(ns("new_name"), "Sample name:", placeholder = "e.g. PatientA_T1"),
      actionButton(
        ns("add_path_btn"),
        "➕ Add Path",
        class = "btn btn-sm btn-secondary",
        style = "width: fit-content;"
      )
    ),
    br(),
    uiOutput(ns("path_list_ui")),
    tags$hr(),
    
    h3("Per-sample metadata"),
    div(
      style = "display: flex; flex-direction: column; gap: 6px; max-width: 500px;",
      textInput(
        ns("new_meta_col"),
        "Add metadata column (letters/numbers/underscore; start with a letter):",
        placeholder = "e.g. sex, condition, batch, donor"
      ),
      actionButton(
        ns("add_meta_col_btn"),
        "➕ Add Column",
        class = "btn btn-sm btn-secondary",
        style = "width: fit-content;"
      )
    ),
    uiOutput(ns("meta_table_ui")),
    br(),
    
    div(
      style = "margin-top: 15px;",
      actionButton(ns("present_btn"), "📦 Create & Present Seurat Object")
    ),
    verbatimTextOutput(ns("seurat_info"))
  )
}



qcModuleServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
   
    # ---------- state ----------
    if (is.null(app_data$seurat_obj)) app_data$seurat_obj <- NULL
    paths_df <- reactiveVal(data.frame(name = character(), path = character(), stringsAsFactors = FALSE))
    meta_cols <- reactiveVal(character())
    dataset  <- reactiveVal(app_data$seurat_obj)
    observe({ app_data$seurat_obj <- dataset() })
   
    # ---------- helpers ----------
    `%||%` <- function(a, b) if (!is.null(a)) a else b
   
    pick_10x_matrix <- function(x) {
      if (inherits(x, c("dgCMatrix","matrix"))) return(x)
      if (is.list(x)) {
        nms <- names(x) %||% rep("", length(x))
        prio <- c("Gene Expression","GeneExpression","RNA","gex","counts","expression")
        for (nm in prio) {
          hit <- which(tolower(nms) == tolower(nm))
          if (length(hit)) return(x[[hit[1]]])
        }
        sizes <- vapply(x, function(m) if (inherits(m, c("dgCMatrix","matrix"))) nrow(m)*ncol(m) else 0, 0)
        if (max(sizes) > 0) return(x[[which.max(sizes)]])
      }
      stop("Could not pick a gene-expression matrix from this 10X input.")
    }
   
    .collect_meta_df <- reactive({
      df <- paths_df()
      cols <- meta_cols()
      if (!nrow(df) || !length(cols)) {
        return(data.frame(row.names = df$name))  # possibly 0 columns, rownames = samples
      }
      samples <- df$name
      out <- lapply(cols, function(col) {
        vapply(samples, function(s) {
          id <- paste0("meta_", col, "_", s)
          val <- input[[id]]
          if (is.null(val) || !nzchar(val)) NA_character__ else as.character(val)
        }, character(1))
      })
      names(out) <- cols
      md <- as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
      rownames(md) <- samples
      md
    })
   
    # ---------- reset ----------
    observeEvent(app_data$reset_qc, {
      paths_df(data.frame(name = character(), path = character(), stringsAsFactors = FALSE))
      meta_cols(character())
      dataset(NULL)
      updateTextInput(session, "new_path", value = "")
      updateTextInput(session, "new_name", value = "")
      updateTextInput(session, "new_meta_col", value = "")
      app_data$reset_qc <- FALSE
    })
   
    # ---------- add path ----------
    observeEvent(input$add_path_btn, {
      p <- trimws(input$new_path)
      n <- trimws(input$new_name)
      if (!nzchar(p) || !file.exists(p)) { showNotification("❌ Invalid path", type="error"); return() }
      if (!nzchar(n)) { showNotification("❌ Please provide a sample name", type="error"); return() }
      df <- paths_df()
      if (n %in% df$name) { showNotification("⚠️ Sample name already exists.", type="warning"); return() }
      if (p %in% df$path) { showNotification("⚠️ Path already added.", type="warning"); return() }
      df <- rbind(df, data.frame(name = n, path = p, stringsAsFactors = FALSE))
      paths_df(df)
      updateTextInput(session, "new_path", value = "")
      updateTextInput(session, "new_name", value = "")
    })
   
    # ---------- add/remove metadata columns ----------
    observeEvent(input$add_meta_col_btn, {
      col <- trimws(input$new_meta_col)
      if (!nzchar(col)) { showNotification("❌ Column name cannot be empty.", type="error"); return() }
      if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", col)) {
        showNotification("❌ Column must start with a letter and contain only letters, digits, or underscores.", type="error")
        return()
      }
      cols <- meta_cols()
      if (col %in% cols) { showNotification("⚠️ Column already exists.", type="warning"); return() }
      meta_cols(c(cols, col))
      updateTextInput(session, "new_meta_col", value = "")
    })
   
    observe({
      lapply(meta_cols(), function(col) {
        observeEvent(input[[paste0("rmcol_", col)]], ignoreInit = TRUE, {
          meta_cols(setdiff(meta_cols(), col))
        })
      })
    })
   
    # ---------- UI renders ----------
    output$path_list_ui <- renderUI({
      df <- paths_df()
      if (!nrow(df)) return(NULL)
      
      tags$ul(
        lapply(seq_len(nrow(df)), function(i) {
          samp <- df$name[i]
          path <- df$path[i]
          
          tags$li(
            div(
              style = "display:flex; align-items:center; gap:6px;",
              tags$span(HTML(paste0("<b>", htmltools::htmlEscape(samp), "</b>: ", htmltools::htmlEscape(path)))),
              actionLink(ns(paste0("rm_path_", samp)), label = HTML("&times;"), title = paste0("Remove ", samp))
            )
          )
        })
      )
    })
    
    observe({
      df <- paths_df()
      lapply(df$name, function(samp) {
        observeEvent(input[[paste0("rm_path_", samp)]], ignoreInit = TRUE, {
          # Remove the sample from paths_df
          df <- paths_df()
          df <- df[df$name != samp, ]
          paths_df(df)
          
          # Remove associated metadata inputs for this sample
          cols <- meta_cols()
          for (col in cols) {
            input_id <- paste0("meta_", col, "_", samp)
            # Reset the input to NULL or empty string
            updateTextInput(session, input_id, value = "")
          }
          
          showNotification(paste0("Removed sample ", samp, " and cleared its metadata."), type = "message")
        })
      })
    })
    
    
   
    output$meta_table_ui <- renderUI({
      df <- paths_df()
      cols <- meta_cols()
     
      if (!nrow(df)) {
        return(helpText("Add sample paths above to enable per-sample metadata."))
      }
      if (!length(cols)) {
        return(helpText("No metadata columns yet. Add one with the field above (e.g., sex, condition, batch, donor)."))
      }
     
      samples <- df$name
      header <- tags$tr(
        tags$th("Sample"),
        lapply(cols, function(col) {
          tags$th(
            div(style = "display:flex; align-items:center; gap:6px;",
                tags$span(col),
                actionLink(ns(paste0("rmcol_", col)), label = HTML("&times;"), title = paste0("Remove ", col))
            )
          )
        })
      )
     
      body <- lapply(samples, function(s) {
        cells <- lapply(cols, function(col) {
          textInput(ns(paste0("meta_", col, "_", s)), label = NULL, value = "", width = "100%")
        })
        do.call(tags$tr, c(list(tags$td(tags$b(s))), lapply(cells, tags$td)))
      })
     
      tags$div(
        style = "overflow-x:auto;",
        tags$table(
          class = "table table-sm",
          style = "width:auto; border-collapse:collapse;",
          tags$thead(header),
          tags$tbody(body)
        )
      )
    })
   
    output$seurat_info <- renderText({
      obj <- dataset()
      if (is.null(obj)) "❌ No object loaded yet." else paste(
        "✅ Seurat object loaded:\n",
        "Class:", class(obj)[1],
        "\nCells:", ncol(obj),
        "\nGenes:", nrow(obj),
        "\n\nProceed to next step"
      )
    })
   
    # ---------- merge + apply metadata ----------
    observeEvent(input$present_btn, {
      df <- paths_df()
      if (!nrow(df)) { showNotification("❌ No paths added", type="error"); return() }
     
      objs <- list()
      for (i in seq_len(nrow(df))) {
        path <- df$path[i]
        samp <- df$name[i]
        
        if (samp %in% names(objs)) {
          showNotification(paste("⚠️ Sample", samp, "already loaded, skipping"), type="warning")
          next
        }
        
        tryCatch({
          raw <- if (grepl("\\.h5$", path, ignore.case = TRUE)) {
            Read10X_h5(path, use.names = TRUE)
          } else {
            Read10X(data.dir = path)
          }
          mat <- pick_10x_matrix(raw)
          # --- Check orientation: genes × cells ---
          if (nrow(mat) < ncol(mat)) {
            warning(sprintf("Matrix for sample %s has (%d rows, %d cols)...", 
                            samp, nrow(mat), ncol(mat)))
          }
          
          # Compute nFeature_RNA and nCount_RNA the Seurat way
          nFeature_RNA <- Matrix::colSums(mat > 0)
          nCount_RNA <- Matrix::colSums(mat)
          
          if (ncol(mat) > 20000) {
            feature_filter <- 200
          } else {
            feature_filter <- 0
          }
          
          keep <- nFeature_RNA > feature_filter 
          mat_filtered <- mat[, keep, drop = FALSE]
          
          
          if (ncol(mat_filtered) == 0) {
            showNotification(paste0("⚠️ All cells filtered out for sample ", samp, ". Skipping this sample."), type="warning")
            next
          }
          
        
          so <- CreateSeuratObject(counts = mat_filtered, project = samp)
        
          so$orig.ident <- samp  # for per-sample metadata
          

          # --- Ensure unique cell names per sample ---
          so <- RenameCells(so, add.cell.id = samp)
          so$orig.ident <- samp
          objs[[samp]] <- so
          
          print(paste("Loaded sample", samp, ":", nrow(so), "genes ×", ncol(so), "cells"))
        }, error = function(e) {
          showNotification(paste0("❌ Error loading ", samp, " from ", path, ": ", e$message),
                           type = "error", duration = 10)
        })
      }
     
      if (!length(objs)) { showNotification("❌ No valid Seurat objects created", type="error"); return() }
     
      for (s in names(objs)) {
        objs[[s]] <- RenameCells(objs[[s]], add.cell.id = s)
      }
      
      
      sample_names <- names(objs)
      merged_obj <- if (length(objs) == 1) {
        so1 <- objs[[1]]
        colnames(so1) <- paste0(sample_names[1], "_", colnames(so1))
        so1
      } else {
        merge(
          x            = objs[[1]],
          y            = objs[-1],
          add.cell.ids = sample_names,
          project      = "Combined"
        )
      }
      
      print(paste("Merged object:", nrow(merged_obj), "genes ×", ncol(merged_obj), "cells across", length(unique(merged_obj$orig.ident)), "samples"))
      
     
      # ---- Apply per-sample metadata using orig.ident ----
      md <- .collect_meta_df()
      cols <- meta_cols()
     
      if (length(cols)) {
        if (nrow(md)) md[md == ""] <- NA
       
        key <- merged_obj$orig.ident
        missing_rows <- setdiff(unique(key), rownames(md))
        if (length(missing_rows)) {
          if (ncol(md) == 0) {
            md <- data.frame(row.names = unique(key))
          } else {
            add <- as.data.frame(
              matrix(NA_character_, nrow = length(missing_rows), ncol = ncol(md),
                     dimnames = list(missing_rows, colnames(md))),
              stringsAsFactors = FALSE
            )
            md <- rbind(md, add)
          }
        }
       
        for (col in cols) {
          vals <- md[key, col]
          suppressWarnings(vnum <- as.numeric(vals))
          if (all(is.na(vals) == is.na(vnum))) {
            merged_obj[[col]] <- vnum
          } else {
            if (length(vals) != ncol(merged_obj)) {
              warning(sprintf(
                "Metadata column '%s' length mismatch: %d vs %d cells. Skipping this column.",
                col, length(vals), ncol(merged_obj)
              ))
              next  # skip this metadata column
            }
            
            merged_obj[[col]] <- vals
          }
        }
      }
     
      showNotification(
        paste0("✅ Merged: ",
               nrow(merged_obj), " genes × ", ncol(merged_obj),
               " cells across ", length(unique(merged_obj$orig.ident)), " samples",
               if (length(cols)) paste0(" | Added metadata: ", paste(cols, collapse = ", ")) else ""),
        type = "message", duration = 8
      )
     
      dataset(merged_obj)
    })
   
    return(dataset)
  })
}
