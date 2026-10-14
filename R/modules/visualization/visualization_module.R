# Cardinal NEST - visualization module
# Ion-image/spectrum viewers plus a Figure Builder that collects displayed
# plots into a grid figure for PNG/PDF download. Metadata image mode renders
# a categorical pixelData column as a spatial map.

library(shiny)
library(ggplot2)

# Panel labels "(a)", "(b)", ..., continuing "(aa)", "(ab)", ... past 26.
collection_panel_labels <- function(n) {
  if (n <= 0) return(character(0))
  base_letters <- letters
  make_label <- function(i) {
    if (i <= 26L) return(base_letters[i])
    # zero-based conversion to a bijective base-26 sequence (aa, ab, ...).
    first <- ((i - 1L) %/% 26L)
    second <- ((i - 1L) %% 26L) + 1L
    paste0(base_letters[first], base_letters[second])
  }
  vapply(seq_len(n), function(i) sprintf("(%s)", make_label(i)), character(1))
}

# Grid dims for n panels; NA ncol gives a near-square layout.
collection_grid_dims <- function(n, ncol = NA) {
  if (n <= 0) return(list(nrow = 0L, ncol = 0L))
  if (is.null(ncol) || is.na(ncol) || ncol <= 0) {
    ncol <- ceiling(sqrt(n))
  }
  ncol <- min(as.integer(ncol), n)
  nrow <- as.integer(ceiling(n / ncol))
  list(nrow = nrow, ncol = ncol)
}

# One single-panel plot per entry in a displayed plot's $plots list.
split_display_plot <- function(p) {
  if (is.null(p)) return(list())
  plots <- p$plots
  n <- if (is.null(plots)) 0L else length(plots)
  if (n <= 1L) return(list(p))
  lapply(seq_len(n), function(i) {
    q <- p
    q$plots <- plots[i]
    if (!is.null(p$labels)) q$labels <- p$labels[i]
    q$subscripts <- NULL
    q$dim <- c(1L, 1L)
    q
  })
}

# Capture one plot as a grob; ggplot via ggplotGrob, Cardinal base-graphics
# via echoGrob (keeps each panel's own coordinates). Optional header label.
collection_panel_grob <- function(plot, label = NULL) {
  if (inherits(plot, "ggplot")) {
    g <- ggplot2::ggplotGrob(plot)
  } else {
    g <- gridGraphics::echoGrob(function() print(plot))
  }
  if (!is.null(label) && nzchar(label)) {
    g <- gridExtra::arrangeGrob(
      g,
      top = grid::textGrob(label, x = grid::unit(0.02, "npc"), hjust = 0,
                           gp = grid::gpar(fontface = "bold", cex = 1.1))
    )
  }
  g
}

# Arrange collected list(plot, label) items into one grob; NULL when empty.
build_plot_collection <- function(items, ncol = NA, labels = TRUE) {
  if (length(items) == 0) return(NULL)
  dims <- collection_grid_dims(length(items), ncol = ncol)
  grobs <- lapply(items, function(it) {
    lab <- if (isTRUE(labels)) it$label else NULL
    collection_panel_grob(it$plot, lab)
  })
  gridExtra::arrangeGrob(grobs = grobs, nrow = dims$nrow, ncol = dims$ncol)
}

# Render the collection to PNG/PDF; pixel dims, converted at 100 dpi for PDF.
render_collection_file <- function(collection, file,
                                    device = c("png", "pdf"),
                                    ncol = NA, labels = TRUE,
                                    width = 1200, height = 900, res = 120) {
  device <- match.arg(device)
  if (is.null(collection)) stop("The plot collection is empty.", call. = FALSE)
  if (is.list(collection) && !inherits(collection, "grob")) {
    if (length(collection) == 0) stop("The plot collection is empty.", call. = FALSE)
    collection <- build_plot_collection(collection, ncol = ncol, labels = labels)
  }
  if (device == "png") {
    grDevices::png(file, width = width, height = height, res = res)
  } else {
    grDevices::pdf(file, width = width / 100, height = height / 100)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  # echoGrob emits harmless is.na() warnings on expression axis labels.
  suppressWarnings(grid::grid.draw(collection))
  invisible(file)
}

# Finite y-range for base graphics `ylim`; NULL when nothing finite.
finite_spectrum_ylim <- function(values) {
  vals <- suppressWarnings(as.numeric(values))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NULL)
  yr <- range(vals)
  if (yr[1] == yr[2]) {
    pad <- abs(yr[1]) * 0.05 + 1
    yr <- yr + c(-pad, pad)
  }
  yr
}

# Caption for a collected plot from mode + selection inputs.
make_plot_caption <- function(mode, selection, detail = "") {
  head <- if (identical(mode, "spectra")) "Spectra" else "Ion image"
  sel <- switch(selection,
                index = "index", coord = "coord", summary = "summary",
                mz = "m/z", selection)
  detail <- trimws(detail)
  if (nzchar(detail)) sprintf("%s \u2014 %s: %s", head, sel, detail)
  else sprintf("%s \u2014 %s", head, sel)
}

# Categorical pixelData columns usable for a spatial map (factors, chars,
# logicals, low-cardinality numerics); coordinate columns excluded.
metadata_candidate_columns <- function(data, max_levels = 50) {
  if (is.null(data)) return(character(0))
  pd <- tryCatch(as.data.frame(Cardinal::pData(data)), error = function(e) NULL)
  if (is.null(pd) || ncol(pd) == 0) return(character(0))
  # Pixel coordinate columns (x/y/z) are spatial, not categorical metadata.
  coord_names <- tryCatch(names(as.data.frame(Cardinal::coord(data))),
                          error = function(e) character(0))
  keep <- vapply(pd, function(col) {
    vals <- col[!is.na(col)]
    nu <- length(unique(vals))
    if (is.factor(col) || is.character(col) || is.logical(col)) {
      return(nu >= 1 && nu <= max_levels)
    }
    if (is.numeric(col)) {
      return(nu >= 1 && nu <= max_levels && nu < length(col))
    }
    FALSE
  }, logical(1))
  setdiff(names(pd)[keep], coord_names)
}

# Spatial map colored by a categorical column; 1-D strip fallback.
plot_metadata_image <- function(coord, values, column = "metadata") {
  values <- as.factor(values)
  n <- nlevels(values)
  palette <- grDevices::hcl.colors(max(n, 2), "Dark 3")
  main <- sprintf("%s (%d categor%s)", column, n, if (n == 1) "y" else "ies")

  if (!is.null(coord) && all(c("x", "y") %in% names(coord)) &&
      nrow(coord) == length(values)) {
    df <- data.frame(coord, category = values)
    ggplot(df, aes(x = x, y = y, fill = category)) +
      geom_tile() +
      coord_equal() +
      scale_fill_manual(values = palette, name = column) +
      labs(title = main, x = "x", y = "y") +
      theme_minimal()
  } else {
    df <- data.frame(pixel = seq_along(values), category = values)
    ggplot(df, aes(x = pixel, y = 1, fill = category)) +
      geom_tile() +
      scale_fill_manual(values = palette, name = column) +
      labs(title = main, x = "pixel", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

visualizationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
             h3("Visualization"),
             p("Plot spectra and render ion images from your MSI dataset"),
             uiOutput(ns("active_dataset"))
      )
    ),
    
    fluidRow(
      column(6,
        bslib::card(
          bslib::card_header("Specifications"),
          bslib::card_body(
            
            selectInput(ns("plot_mode"), "Plot Mode",
                        choices = c("Ion Image" = "ion",
                                    "Spectra" = "spectra",
                                    "Metadata Image" = "metadata"),
                        selected = "spectra"),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'metadata'", ns("plot_mode")),
              selectInput(ns("meta_column"), "Categorical column",
                          choices = character(0)),
              p(class = "text-muted small",
                "Color pixels by a categorical column such as group labels, ",
                "clusters, segments, or a saved ROI.")
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'spectra'", ns("plot_mode")),
              
              selectInput(ns("spec_select"), "Spectrum Selection",
                          choices = c("Index" = "index",
                                      "Pixel Coordinate" = "coord",
                                      "Summary" = "summary")),
              
              selectInput(ns("spec_style"), "Spectrum Style",
                          choices = c("Histogram" = "h",
                                      "Line" = "l"),
                          selected = "h"),
              
              conditionalPanel(
                sprintf("input['%s'] == 'index'", ns("spec_select")),
                textInput(ns("spec_index"), "Spectrum Index", value = "",
                          placeholder = "e.g. 1:10, 30")
              ),
              
              conditionalPanel(
                sprintf("input['%s'] == 'coord'", ns("spec_select")),
                numericInput(ns("spec_x"), "Pixel X Coordinate", value = NA),
                numericInput(ns("spec_y"), "Pixel Y Coordinate", value = NA)
              ),
              
              conditionalPanel(
                sprintf("input['%s'] == 'summary'", ns("spec_select")),
                selectInput(ns("spec_summary"), "Summary Selection",
                            choices = c("Mean" = "mean"))
              ),
              
              checkboxInput(ns("spec_superpose"), "Superpose", value = FALSE),
              
              shinyWidgets::numericRangeInput(ns("spec_mz_range"), "m/z range",
                                              value = c(NA, NA), separator = " to "),
              
              bslib::accordion(
                open = FALSE,
                bslib::accordion_panel(
                  "Advanced features",
                  numericInput(ns("spec_line_width"), "Line Width", value = NA, step = "any"),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'summary'", ns("spec_select")),
                    numericInput(ns("spec_ann_peaks"), "Annotations", value = NA)
                  ),
                  checkboxInput(ns("spec_show_key"), "Key", value = FALSE),
                  textInput(ns("spec_xlab"), "X Label", value = "m/z", width = "100%",
                            placeholder = "e.g., m/z"),
                  textInput(ns("spec_ylab"), "Y Label", value = "Intensity", width = "100%",
                            placeholder = "e.g., Intensity")
                )
              )
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'ion'", ns("plot_mode")),
              
              selectInput(ns("ion_select"), "Ion Selection",
                          choices = c("m/z" = "mz",
                                      "Feature Index" = "index",
                                      "Summary" = "summary")),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'mz'", ns("ion_select")),
                textInput(ns("ion_mz"), "m/z", value = "",
                          placeholder = "e.g. 1003.3 or 1003.3, 1663.6")
              ),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'index'", ns("ion_select")),
                textInput(ns("ion_i"), "Feature Index", value = "",
                          placeholder = "e.g. 1:10, 30")
              ),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'summary'", ns("ion_select")),
                selectInput(ns("ion_summary_select"), "Pixel Summary",
                            choices = c("TIC" = "sum",
                                        "Mean" = "mean"))
              ),
              
              checkboxInput(ns("ion_superpose"), "Superpose", value = FALSE),
              
              bslib::accordion(
                open = FALSE,
                bslib::accordion_panel(
                  "Advanced features",
                  numericInput(ns("ion_tolerance"), "Tolerance", value = 0.5, step = "any"),
                  selectInput(ns("ion_units"), "Units",
                              choices = c("m/z" = "mz", "ppm" = "ppm"), selected = "mz"),
                  selectInput(ns("ion_smooth"), "Smooth",
                              choices = c("None" = "", "Gaussian" = "gaussian",
                                          "Adaptive" = "adaptive", "Guide" = "guide")),
                  selectInput(ns("ion_enhance"), "Enhance",
                              choices = c("None" = "", "Histogram" = "histogram",
                                          "Adaptive" = "adaptive")),
                  checkboxInput(ns("ion_scale"), "Scale", value = FALSE),
                  textInput(ns("ion_run"), "Run", value = "", placeholder = "e.g. runA1"),
                  selectInput(ns("ion_color"), "Color Palette",
                              choices = c("Default" = "", "Cividis" = "Cividis",
                                          "Viridis" = "viridis", "Plasma" = "plasma")),
                  selectInput(ns("ion_style"), "Style",
                              choices = c("Default" = "", "Dark" = "dark", "Light" = "light"))
                )
              )
            ),
            
            div(class = "d-flex gap-2",
              actionButton(ns("display_image"), "Display", class = "btn-primary"),
              actionButton(ns("add_to_collection"), "Add to collection",
                           icon = icon("plus"), class = "btn-outline-success")
            )
          )
        )
      ),
      
      column(6,
        bslib::card(
          bslib::card_header("Result"),
          bslib::card_body(
            shinycssloaders::withSpinner(
              plotOutput(ns("viz_result"), height = "75vh")
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header(
            div(class = "d-flex justify-content-between align-items-center",
              span(icon("table-cells"), " Figure Builder"),
              uiOutput(ns("collection_count"))
            )
          ),
          bslib::card_body(
            p(class = "text-muted small",
              "Display a spectrum or ion image above, then click ",
              strong("Add to collection"), " to place it in the figure grid. ",
              "Plots shown as multiple panels are added as separate, individually ",
              "labelled items. Use ", strong("Preview"),
              " to see the assembled figure before downloading."),
            fluidRow(
              column(4,
                numericInput(ns("coll_ncol"), "Columns (blank = auto)",
                             value = NA, min = 1, step = 1)
              ),
              column(4,
                div(style = "margin-top: 2rem;",
                  checkboxInput(ns("coll_labels"),
                                "Show panel labels", value = TRUE)
                )
              )
            ),
            div(class = "d-flex gap-2 mb-3",
              actionButton(ns("preview_collection"), "Preview",
                           icon = icon("eye"),
                           class = "btn-outline-secondary btn-sm"),
              downloadButton(ns("download_png"), "PNG",
                             class = "btn-outline-primary btn-sm"),
              downloadButton(ns("download_pdf"), "PDF",
                             class = "btn-outline-primary btn-sm"),
              actionButton(ns("clear_collection"), "Clear all",
                           icon = icon("trash"),
                           class = "btn-outline-danger btn-sm")
            ),
            uiOutput(ns("collection_items"))
          )
        )
      )
    )
  )
}

visualizationServer <- function(id, datasets, active_id) {
  moduleServer(id, function(input, output, session) {
    
    viz_result <- reactiveVal(NULL)
    error_msg <- reactiveVal(NULL)
    current_caption <- reactiveVal(NULL)

    # Ordered collected plots: each list(uid, plot, caption, label).
    plot_collection <- reactiveVal(list())
    collection_counter <- reactiveVal(0)
    
    mse_data <- active_dataset(datasets, active_id)

    # parse_numeric_list wrapper: NULL + recorded error on invalid input.
    parse_numeric_input <- function(str, label, integer = FALSE) {
      tryCatch(parse_numeric_list(str, label, integer = integer),
               error = function(e) {
                 viz_result(NULL)
                 error_msg(conditionMessage(e))
                 NULL
               })
    }
    
    # Shared progress + error handling around a plotting call.
    render_plot <- function(fn, args, post = identity) {
      tryCatch({
        p <- withProgress(message = "Rendering plot", value = 0.2, {
          setProgress(0.5, detail = "Computing")
          out <- post(do.call(fn, args))
          setProgress(1, detail = "Done")
          out
        })
        viz_result(p)
        error_msg(NULL)
      }, error = function(e) {
        viz_result(NULL)
        error_msg(conditionMessage(e))
      })
    }
    
    # Expand spectrum axis limits to "pretty" bounds.
    expand_spectra_axes <- function(p) {
      # pretty() on NA/Inf yields non-finite limits and crashes base graphics.
      expand <- function(lims) {
        finite <- lims[is.finite(lims)]
        if (length(finite) == 0) return(lims)
        range(pretty(range(finite)))
      }
      if (!is.null(p$channels$x$limits))
        p$channels$x$limits <- expand(p$channels$x$limits)
      if (!is.null(p$channels$y$limits))
        p$channels$y$limits <- expand(p$channels$y$limits)
      p
    }
    
    observe({
      if (is.null(mse_data())) {
        shinyjs::disable("display_image")
      } else {
        shinyjs::enable("display_image")
      }
    })

    # Populate the categorical-column picker when the dataset changes.
    observe({
      cols <- metadata_candidate_columns(mse_data())
      updateSelectInput(session, "meta_column",
                        choices = cols,
                        selected = if (length(cols)) cols[[1]] else character(0))
    })
    
    output$active_dataset <- active_dataset_ui(datasets, active_id)

    observeEvent(input$display_image, {
      req(mse_data())
      
      if (input$plot_mode == "spectra") {
        args <- list(mse_data())
        
        if (input$spec_select == "index") {
          idx <- parse_numeric_input(input$spec_index, "Spectrum index input",
                                     integer = TRUE)
          if (is.null(idx)) return()
          args$i <- idx
        } else if (input$spec_select == "coord") {
          req(!is.na(input$spec_x), !is.na(input$spec_y))
          args$coord <- list(x = input$spec_x, y = input$spec_y)
        } else if (input$spec_select == "summary") {
          req(input$spec_summary)
          data_sum <- tryCatch(
            summarizeFeatures(mse_data(), input$spec_summary),
            error = function(e) {
              error_msg(paste("Could not summarize features:", conditionMessage(e)))
              viz_result(NULL)
              NULL
            }
          )
          if (is.null(data_sum)) return()
          # Pin a finite ylim; non-finite summary values crash base graphics.
          svals <- tryCatch(Cardinal::fData(data_sum)[[input$spec_summary]],
                            error = function(e) NULL)
          args <- list(data_sum, input$spec_summary)
          if (length(svals) > 0) {
            summary_ylim <- finite_spectrum_ylim(svals)
            if (is.null(summary_ylim)) {
              viz_result(NULL)
              error_msg(sprintf(
                "The %s summary spectrum has no finite intensity values to plot.",
                input$spec_summary))
              return()
            }
            args$ylim <- summary_ylim
          }
        }
        
        mzr <- input$spec_mz_range
        if (length(mzr) == 2 && all(!is.na(mzr))) {
          args$xlim <- c(mzr[1], mzr[2])
        }
        
        args$superpose <- input$spec_superpose

        if (isTruthy(input$spec_line_width) && identical(input$spec_style, "l"))
          args$linewidth <- input$spec_line_width
        if (input$spec_select == "summary" && isTruthy(input$spec_ann_peaks))
          args$annPeaks <- input$spec_ann_peaks
        args$key <- input$spec_show_key
        args$xlab <- input$spec_xlab
        args$ylab <- input$spec_ylab
        args$type <- if (isTruthy(input$spec_style)) input$spec_style else "h"

        spec_detail <- switch(input$spec_select,
          index = input$spec_index,
          coord = sprintf("(%s, %s)", input$spec_x, input$spec_y),
          summary = input$spec_summary, "")
        current_caption(make_plot_caption("spectra", input$spec_select, spec_detail))

        render_plot(plot, args, post = expand_spectra_axes)

      } else if (input$plot_mode == "metadata") {
        data <- mse_data()
        col <- input$meta_column
        if (!isTruthy(col)) {
          viz_result(NULL)
          error_msg("Select a categorical column to display.")
          return()
        }
        pd <- as.data.frame(Cardinal::pData(data))
        if (!col %in% names(pd)) {
          viz_result(NULL)
          error_msg(sprintf("Column '%s' is not present in the active dataset.", col))
          return()
        }
        coord <- as.data.frame(Cardinal::coord(data))
        current_caption(sprintf("Metadata image \u2014 %s", col))
        render_plot(plot_metadata_image,
                    list(coord = coord, values = pd[[col]], column = col))

      } else {
        args <- list(mse_data())
        
        if (input$ion_select == "mz") {
          mz_vals <- parse_numeric_input(input$ion_mz, "m/z input")
          if (is.null(mz_vals)) return()
          args$mz <- mz_vals
        } else if (input$ion_select == "index") {
          i_vals <- parse_numeric_input(input$ion_i, "Feature index input",
                                        integer = TRUE)
          if (is.null(i_vals)) return()
          args$i <- i_vals
        } else if (input$ion_select == "summary") {
          req(input$ion_summary_select)
          data_sum <- tryCatch(
            summarizePixels(mse_data(), input$ion_summary_select),
            error = function(e) {
              error_msg(paste("Could not summarize pixels:", conditionMessage(e)))
              viz_result(NULL)
              NULL
            }
          )
          if (is.null(data_sum)) return()
          args <- list(data_sum, input$ion_summary_select)
        }
        
        args$superpose <- input$ion_superpose

        if (isTruthy(input$ion_tolerance)) args$tolerance <- input$ion_tolerance
        if (input$ion_units != "") args$units <- input$ion_units
        if (input$ion_smooth != "") args$smooth <- input$ion_smooth
        if (input$ion_enhance != "") args$enhance <- input$ion_enhance
        if (input$ion_scale) args$scale <- input$ion_scale
        if (input$ion_run != "") args$run <- input$ion_run
        if (input$ion_color != "") args$col <- cpal(input$ion_color)
        if (input$ion_style != "") args$style <- input$ion_style

        ion_detail <- switch(input$ion_select,
          mz = input$ion_mz,
          index = input$ion_i,
          summary = input$ion_summary_select, "")
        current_caption(make_plot_caption("ion", input$ion_select, ion_detail))

        render_plot(image, args)
      }
    })
    
    output$viz_result <- renderPlot({
      validate(need(is.null(error_msg()), error_msg()))
      p <- viz_result()
      req(p)
      print(p)
    })

    observe({
      if (is.null(viz_result())) {
        shinyjs::disable("add_to_collection")
      } else {
        shinyjs::enable("add_to_collection")
      }
    })

    observe({
      empty <- length(plot_collection()) == 0
      shinyjs::toggleState("clear_collection", condition = !empty)
      shinyjs::toggleState("preview_collection", condition = !empty)
      shinyjs::toggleState("download_png", condition = !empty)
      shinyjs::toggleState("download_pdf", condition = !empty)
    })

    # Sync live label inputs back into the collection before re-render.
    sync_labels <- function(current = plot_collection()) {
      if (length(current) == 0) return(current)
      for (i in seq_along(current)) {
        val <- input[[paste0("label_", current[[i]]$uid)]]
        if (!is.null(val)) current[[i]]$label <- val
      }
      current
    }

    # Multi-panel plots (non-superposed spectra) split into separate items.
    observeEvent(input$add_to_collection, {
      p <- viz_result()
      if (is.null(p)) return()
      parts <- split_display_plot(p)
      base_cap <- current_caption()
      current <- sync_labels()
      multi <- length(parts) > 1
      for (j in seq_along(parts)) {
        new_count <- collection_counter() + 1
        collection_counter(new_count)
        pos <- length(current) + 1
        def_lab <- collection_panel_labels(pos)[pos]
        cap <- base_cap
        if (is.null(cap)) cap <- sprintf("Plot %d", new_count)
        if (multi) cap <- sprintf("%s [panel %d]", cap, j)
        current[[pos]] <- list(
          uid = paste0("c", new_count), plot = parts[[j]],
          caption = cap, label = def_lab
        )
      }
      plot_collection(current)
    })

    observeEvent(input$remove_plot, {
      uid <- input$remove_plot$uid
      current <- sync_labels()
      keep <- vapply(current, function(it) it$uid != uid, logical(1))
      plot_collection(current[keep])
    })

    observeEvent(input$clear_collection, {
      plot_collection(list())
    })

    # Stored plots paired with live label text (blank when labels are off).
    resolved_items <- reactive({
      items <- plot_collection()
      show <- isTRUE(input$coll_labels)
      lapply(items, function(it) {
        lab <- input[[paste0("label_", it$uid)]]
        if (is.null(lab)) lab <- it$label
        list(plot = it$plot, label = if (show) lab else "")
      })
    })

    # Assembled figure grob shared by the preview and both downloads.
    collection_figure <- reactive({
      items <- resolved_items()
      if (length(items) == 0) return(NULL)
      build_plot_collection(items, ncol = input$coll_ncol, labels = TRUE)
    })

    output$collection_count <- renderUI({
      n <- length(plot_collection())
      span(class = "badge bg-secondary",
           sprintf("%d panel%s", n, if (n == 1) "" else "s"))
    })

    output$collection_items <- renderUI({
      items <- plot_collection()
      if (length(items) == 0) {
        return(div(class = "text-muted small",
                   "No panels yet. Displayed plots you add will appear here."))
      }
      tagList(
        lapply(seq_along(items), function(i) {
          it <- items[[i]]
          div(class = "d-flex align-items-center gap-2 border rounded px-2 py-1 mb-1",
            div(style = "width: 90px;",
              textInput(session$ns(paste0("label_", it$uid)), label = NULL,
                        value = it$label, width = "90px")
            ),
            span(class = "flex-grow-1 small text-muted", it$caption),
            tags$button(
              type = "button",
              class = "btn btn-sm btn-outline-danger",
              title = "Remove",
              onclick = sprintf(
                "Shiny.setInputValue('%s', {uid: '%s', nonce: Math.random()}, {priority: 'event'})",
                session$ns("remove_plot"), it$uid
              ),
              icon("times")
            )
          )
        })
      )
    })

    # Modal preview of the assembled figure (no inline preview).
    observeEvent(input$preview_collection, {
      if (length(plot_collection()) == 0) return()
      showModal(modalDialog(
        title = "Figure preview",
        size = "xl", easyClose = TRUE,
        shinycssloaders::withSpinner(
          plotOutput(session$ns("preview_modal_plot"), height = "70vh")
        ),
        footer = modalButton("Close")
      ))
    })

    output$preview_modal_plot <- renderPlot({
      fig <- collection_figure()
      validate(need(!is.null(fig), "Add one or more plots to build a figure."))
      grid::grid.newpage()
      # Harmless is.na() warnings while echoGrob replays expression labels.
      suppressWarnings(grid::grid.draw(fig))
    })

    # Grid dimensions scale the output size.
    collection_dims <- reactive({
      collection_grid_dims(length(plot_collection()), ncol = input$coll_ncol)
    })

    output$download_png <- downloadHandler(
      filename = function() paste0("figure_", Sys.Date(), ".png"),
      content = function(file) {
        fig <- collection_figure()
        req(fig)
        dims <- collection_dims()
        render_collection_file(fig, file, device = "png",
                               width = 500 * dims$ncol,
                               height = 400 * dims$nrow)
      }
    )

    output$download_pdf <- downloadHandler(
      filename = function() paste0("figure_", Sys.Date(), ".pdf"),
      content = function(file) {
        fig <- collection_figure()
        req(fig)
        dims <- collection_dims()
        render_collection_file(fig, file, device = "pdf",
                               width = 500 * dims$ncol,
                               height = 400 * dims$nrow)
      }
    )

  })
}
