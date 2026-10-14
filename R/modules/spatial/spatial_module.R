# Cardinal NEST - spatial module
# SSC/spatial-DGMM segmentation and ROI selection; results can be saved back
# as new pixelData columns on the dataset.

library(shiny)

# Cardinal::spatialShrunkenCentroids wrapper. r = neighbourhood radius,
# k = initial segment count, s = shrinkage, weights = "gaussian"|"adaptive".
# `segment` is per-pixel; `statistic` is feature x segment t-statistics.
run_spatial_segmentation <- function(data, r = 1, k = 3, s = 0,
                                      weights = c("gaussian", "adaptive")) {
  if (is.null(data)) stop("No dataset supplied for segmentation.", call. = FALSE)
  weights <- match.arg(weights)

  k <- as.integer(k)
  if (is.na(k) || k < 2) k <- 2L
  if (ncol(data) < 2L) {
    stop("Dataset has too few pixels for segmentation (need at least 2).",
         call. = FALSE)
  }
  k <- min(k, ncol(data))

  if (is.null(r) || is.na(r) || r < 0) r <- 1
  if (is.null(s) || is.na(s) || s < 0) s <- 0
  r <- as.numeric(r)
  s <- as.numeric(s)

  model <- Cardinal::spatialShrunkenCentroids(data, r = r, k = k, s = s,
                                              weights = weights)

  segment <- as.factor(model$class)
  segment <- factor(segment, levels = sort(unique(as.character(segment))))

  list(
    model = model,
    r = r,
    k = k,
    s = s,
    weights = weights,
    segment = segment,
    n_segments = nlevels(segment),
    sizes = table(segment),
    coord = as.data.frame(Cardinal::coord(data)),
    statistic = model$statistic,
    mz = Cardinal::mz(data)
  )
}

# Top distinguishing features per segment via Cardinal::topFeatures().
segmentation_top_features <- function(seg, n = 20) {
  if (is.null(seg$model)) {
    return(data.frame(i = integer(0), mz = numeric(0),
                      class = character(0), statistic = numeric(0),
                      centers = numeric(0), sd = numeric(0)))
  }
  as.data.frame(Cardinal::topFeatures(seg$model, n = n, sort.by = "statistic"))
}

# Per-pixel segment + coordinates (CSV export / save-back base table).
segmentation_pixel_table <- function(seg) {
  cbind(seg$coord, Segment = as.integer(seg$segment))
}

# Segment map by pixel coordinates; 1-D strip fallback.
plot_segmentation <- function(seg) {
  coord <- seg$coord
  segs <- seg$segment
  palette <- grDevices::hcl.colors(max(nlevels(segs), 2), "Dark 3")
  main <- sprintf("Spatial segmentation (r = %g, k = %d, s = %g) -> %d segments",
                  seg$r, seg$k, seg$s, nlevels(segs))

  if (all(c("x", "y") %in% names(coord))) {
    plot(coord$x, coord$y, pch = 15, col = palette[as.integer(segs)],
         xlab = "x", ylab = "y", main = main, asp = 1)
  } else {
    plot(seq_along(segs), rep(1, length(segs)), pch = 15,
         col = palette[as.integer(segs)],
         xlab = "pixel", ylab = "", yaxt = "n", main = main)
  }
  legend("topright", legend = paste("Segment", levels(segs)),
         col = palette, pch = 15, bty = "n")
  invisible(TRUE)
}

# Dataset record with the segments as pixelData column spatial_segment.
make_segmentation_record <- function(data, seg, id, filename = NULL, ext = NULL) {
  make_result_record(data, list(spatial_segment = seg$segment), id,
                     sprintf("Spatial segmentation (%d segments)", seg$n_segments),
                     filename = filename, ext = ext)
}

# Per-pixel class for `feature`; model$class may be a list, matrix, or
# single factor.
dgmm_feature_class <- function(model, feature = 1) {
  cls <- model$class
  feature <- as.integer(feature)
  if (is.na(feature) || feature < 1) feature <- 1L
  if (is.list(cls) && !is.data.frame(cls)) {
    feature <- min(feature, length(cls))
    v <- cls[[feature]]
  } else if (is.matrix(cls) || is.data.frame(cls)) {
    feature <- min(feature, ncol(cls))
    v <- cls[, feature]
  } else {
    v <- cls
  }
  as.factor(v)
}

# Cardinal::spatialDGMM wrapper; DGMM fits one mixture per feature so a
# single 1-based `feature` is chosen for the per-pixel segment map.
run_spatial_dgmm <- function(data, r = 1, k = 3, feature = 1,
                             weights = c("gaussian", "adaptive")) {
  if (is.null(data)) stop("No dataset supplied for spatial-DGMM.", call. = FALSE)
  weights <- match.arg(weights)

  k <- as.integer(k)
  if (is.na(k) || k < 2) k <- 2L
  if (ncol(data) < 2L) {
    stop("Dataset has too few pixels for spatial-DGMM (need at least 2).",
         call. = FALSE)
  }

  if (is.null(r) || is.na(r) || r < 0) r <- 1
  r <- as.numeric(r)

  n_feat <- nrow(data)
  feature <- as.integer(feature)
  if (is.na(feature) || feature < 1) feature <- 1L
  feature <- min(feature, n_feat)

  model <- Cardinal::spatialDGMM(data, r = r, k = k, weights = weights)

  segment <- dgmm_feature_class(model, feature)
  keep <- as.character(segment[!is.na(segment)])
  segment <- factor(as.character(segment), levels = sort(unique(keep)))

  mz <- Cardinal::mz(data)
  mz_val <- if (!is.null(mz) && length(mz) >= feature) mz[feature] else NA_real_

  list(
    model = model,
    r = r,
    k = k,
    weights = weights,
    feature = feature,
    mz = mz_val,
    segment = segment,
    n_segments = nlevels(segment),
    sizes = table(segment),
    coord = as.data.frame(Cardinal::coord(data))
  )
}

# Per-pixel DGMM segment + coordinates (CSV export / save-back base table).
dgmm_pixel_table <- function(dgmm) {
  cbind(dgmm$coord, Segment = as.character(dgmm$segment))
}

# Cardinal's own image(model, i = feature) rendering of the segment map.
plot_dgmm <- function(dgmm) {
  print(Cardinal::image(dgmm$model, i = dgmm$feature))
  invisible(TRUE)
}

# Dataset record with the DGMM segments as pixelData column dgmm_segment.
make_dgmm_record <- function(data, dgmm, id, filename = NULL, ext = NULL) {
  mz_lab <- if (is.na(dgmm$mz)) sprintf("feature %d", dgmm$feature) else sprintf("m/z %.4g", dgmm$mz)
  make_result_record(data, list(dgmm_segment = dgmm$segment), id,
                     sprintf("Spatial-DGMM (%s, %d segments)", mz_lab, dgmm$n_segments),
                     filename = filename, ext = ext)
}

# Logical vector: pixels inside the brush rectangle.
pixels_in_brush <- function(coord, brush) {
  n <- if (is.null(coord)) 0L else nrow(coord)
  if (n == 0L) return(logical(0))
  if (is.null(brush) || !all(c("x", "y") %in% names(coord))) return(logical(n))
  if (is.null(brush$xmin) || is.null(brush$xmax) ||
      is.null(brush$ymin) || is.null(brush$ymax)) return(logical(n))
  coord$x >= brush$xmin & coord$x <= brush$xmax &
    coord$y >= brush$ymin & coord$y <= brush$ymax
}

# Logical vector: pixels inside polygon `poly` (even-odd ray casting).
pixels_in_polygon <- function(coord, poly) {
  n <- if (is.null(coord)) 0L else nrow(coord)
  if (n == 0L) return(logical(0))
  if (is.null(poly) || !all(c("x", "y") %in% names(coord)) ||
      !all(c("x", "y") %in% names(poly))) return(logical(n))
  px <- poly$x; py <- poly$y
  nv <- length(px)
  if (nv < 3L) return(logical(n))
  x <- coord$x; y <- coord$y
  inside <- rep(FALSE, n)
  j <- nv
  for (i in seq_len(nv)) {
    crosses <- ((py[i] > y) != (py[j] > y)) &
      (x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i])
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

# Per-pixel ROI factor: roi_name for selected, "background" otherwise.
make_roi_labels <- function(coord, selected, roi_name = "roi") {
  n <- if (is.null(coord)) 0L else nrow(coord)
  if (is.null(selected) || length(selected) != n) selected <- rep(FALSE, n)
  roi_name <- if (is.null(roi_name) || !nzchar(roi_name)) "roi" else roi_name
  bg <- if (identical(roi_name, "background")) "not_roi" else "background"
  labels <- ifelse(selected, roi_name, bg)
  factor(labels, levels = c(roi_name, bg))
}

# ROI membership paired with coordinates for CSV export.
roi_pixel_table <- function(coord, selected, roi_name = "roi") {
  labels <- make_roi_labels(coord, selected, roi_name)
  cbind(coord, ROI = as.character(labels))
}

# Pixel map for ROI selection; selected pixels highlighted, optional
# in-progress `poly` outline overlaid (dashed close once 3+ vertices).
plot_roi <- function(coord, selected = NULL, roi_name = "roi", poly = NULL) {
  n <- if (is.null(coord)) 0L else nrow(coord)
  if (is.null(selected) || length(selected) != n) selected <- rep(FALSE, n)
  sel_col <- "#D7191C"
  bg_col <- "#BDBDBD"
  poly_col <- "#2C7FB8"
  cols <- ifelse(selected, sel_col, bg_col)
  main <- sprintf("ROI selection: %d of %d pixels", sum(selected), n)

  has_xy <- all(c("x", "y") %in% names(coord))
  if (has_xy) {
    plot(coord$x, coord$y, pch = 15, col = cols,
         xlab = "x", ylab = "y", main = main, asp = 1)
  } else {
    plot(seq_len(n), rep(1, n), pch = 15, col = cols,
         xlab = "pixel", ylab = "", yaxt = "n", main = main)
  }

  # Overlay the in-progress polygon (x/y maps only).
  if (has_xy && !is.null(poly) && all(c("x", "y") %in% names(poly)) &&
      length(poly$x) >= 1L) {
    graphics::points(poly$x, poly$y, pch = 19, col = poly_col, cex = 1.2)
    if (length(poly$x) >= 2L)
      graphics::lines(poly$x, poly$y, col = poly_col, lwd = 2)
    if (length(poly$x) >= 3L)
      graphics::lines(poly$x[c(length(poly$x), 1L)],
                      poly$y[c(length(poly$y), 1L)],
                      col = poly_col, lwd = 2, lty = 2)
  }

  legend("topright", legend = c(roi_name, "background"),
         col = c(sel_col, bg_col), pch = 15, bty = "n")
  invisible(TRUE)
}

# Dataset record with ROI membership as a pixelData column (default "roi").
make_roi_record <- function(data, selected, id, roi_name = "roi",
                            column = NULL, filename = NULL, ext = NULL) {
  coord <- as.data.frame(Cardinal::coord(data))
  labels <- make_roi_labels(coord, selected, roi_name)
  column <- if (is.null(column) || !nzchar(column)) "roi" else column
  make_result_record(data, stats::setNames(list(labels), column), id,
                     sprintf("ROI '%s' (%d pixels)", roi_name, sum(labels == roi_name)),
                     filename = filename, ext = ext)
}

spatialUI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
             h3("Spatial Analysis"),
             p("Segment the active dataset into spatially-coherent regions or select a region of interest (ROI) directly on the pixel map."),
             uiOutput(ns("active_dataset"))
      )
    ),

    bslib::navset_tab(
      bslib::nav_panel("Segmentation",
        fluidRow(
          column(4,
            bslib::card(
              bslib::card_header("Spatial segmentation"),
              bslib::card_body(
                div(class = "text-muted small mb-2",
                    "Spatial Shrunken Centroids (SSC) partitions pixels into ",
                    "spatially-smoothed segments and ranks the features that ",
                    "distinguish them."),

                numericInput(ns("seg_k"), "Number of segments (k)",
                             value = 3, min = 2, step = 1),
                numericInput(ns("seg_r"), "Spatial neighborhood radius (r)",
                             value = 1, min = 0, step = 1),
                numericInput(ns("seg_s"), "Shrinkage (s)",
                             value = 0, min = 0, step = 1),
                div(class = "text-muted small mb-2",
                    "Larger shrinkage removes uninformative features and tends to ",
                    "merge similar segments."),
                selectInput(ns("seg_weights"), "Spatial weights",
                            choices = c("Gaussian" = "gaussian",
                                        "Adaptive" = "adaptive")),

                actionButton(ns("run_segmentation"), "Run segmentation",
                             icon = icon("play"), class = "btn-success w-100")
              )
            )
          ),

          column(8,
            bslib::card(
              bslib::card_header("Status"),
              bslib::card_body(uiOutput(ns("status")))
            ),
            bslib::card(
              bslib::card_header("Results"),
              bslib::card_body(
                bslib::navset_tab(
                  bslib::nav_panel("Segment map",
                    shinycssloaders::withSpinner(plotOutput(ns("result_plot"), height = "60vh"))
                  ),
                  bslib::nav_panel("Top features",
                    DT::DTOutput(ns("result_table"))
                  )
                )
              )
            ),
            bslib::card(
              bslib::card_header("Export"),
              bslib::card_body(
                fluidRow(
                  column(4, downloadButton(ns("export_csv"), "Download segments (.csv)",
                                           class = "btn-outline-primary w-100")),
                  column(4,
                    selectInput(ns("plot_format"), NULL,
                                choices = c("PNG" = "png", "PDF" = "pdf"), width = "100%"),
                    downloadButton(ns("export_plot"), "Download plot",
                                   class = "btn-outline-primary w-100")
                  ),
                  column(4, uiOutput(ns("save_segments_ui")))
                )
              )
            )
          )
        )
      ),

      bslib::nav_panel("DGMM segmentation",
        fluidRow(
          column(4,
            bslib::card(
              bslib::card_header("Spatial-DGMM segmentation"),
              bslib::card_body(
                div(class = "text-muted small mb-2",
                    "Spatial Dirichlet Gaussian Mixture modelling (spatial-DGMM) ",
                    "fits a probabilistic segmentation to each mass feature ",
                    "individually. Choose a feature (m/z) to view its segment map. ",
                    "This complements SSC, which segments using all features together."),

                numericInput(ns("dgmm_k"), "Gaussian components (k)",
                             value = 3, min = 2, step = 1),
                numericInput(ns("dgmm_r"), "Spatial neighborhood radius (r)",
                             value = 1, min = 0, step = 1),
                numericInput(ns("dgmm_mz"), "Feature m/z (blank = first feature)",
                             value = NA, step = "any"),
                selectInput(ns("dgmm_weights"), "Spatial weights",
                            choices = c("Gaussian" = "gaussian",
                                        "Adaptive" = "adaptive")),
                div(class = "text-muted small mb-2",
                    "spatial-DGMM fits one model per feature and can be slow on ",
                    "large datasets."),

                actionButton(ns("run_dgmm"), "Run spatial-DGMM",
                             icon = icon("play"), class = "btn-success w-100")
              )
            )
          ),

          column(8,
            bslib::card(
              bslib::card_header("Status"),
              bslib::card_body(uiOutput(ns("dgmm_status")))
            ),
            bslib::card(
              bslib::card_header("Results"),
              bslib::card_body(
                bslib::navset_tab(
                  bslib::nav_panel("Segment map",
                    shinycssloaders::withSpinner(plotOutput(ns("dgmm_plot"), height = "60vh"))
                  ),
                  bslib::nav_panel("Pixel segments",
                    DT::DTOutput(ns("dgmm_table"))
                  )
                )
              )
            ),
            bslib::card(
              bslib::card_header("Export"),
              bslib::card_body(
                fluidRow(
                  column(4, downloadButton(ns("dgmm_export_csv"), "Download segments (.csv)",
                                           class = "btn-outline-primary w-100")),
                  column(4,
                    selectInput(ns("dgmm_plot_format"), NULL,
                                choices = c("PNG" = "png", "PDF" = "pdf"), width = "100%"),
                    downloadButton(ns("dgmm_export_plot"), "Download plot",
                                   class = "btn-outline-primary w-100")
                  ),
                  column(4, uiOutput(ns("save_dgmm_ui")))
                )
              )
            )
          )
        )
      ),

      bslib::nav_panel("ROI selection",
        fluidRow(
          column(4,
            bslib::card(
              bslib::card_header("Region of interest"),
              bslib::card_body(
                radioButtons(ns("roi_sel_mode"), "Selection shape",
                             choices = c("Rectangle" = "rect",
                                         "Freehand polygon" = "poly"),
                             selected = "rect", inline = TRUE),

                conditionalPanel(
                  condition = sprintf("input['%s'] == 'rect'", ns("roi_sel_mode")),
                  div(class = "text-muted small mb-2",
                      "Drag a rectangle on the map to select pixels, then add or ",
                      "remove them from the ROI. Repeat to build up complex regions.")
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'poly'", ns("roi_sel_mode")),
                  div(class = "text-muted small mb-2",
                      "Click points on the map to trace an irregular outline ",
                      "(3+ vertices), then add or remove the enclosed pixels."),
                  div(class = "d-flex gap-2 mb-2",
                    actionButton(ns("roi_poly_undo"), "Undo point",
                                 icon = icon("rotate-left"),
                                 class = "btn-outline-secondary btn-sm w-100"),
                    actionButton(ns("roi_poly_clear"), "Clear outline",
                                 icon = icon("eraser"),
                                 class = "btn-outline-secondary btn-sm w-100")
                  ),
                  div(class = "text-muted small mb-2", uiOutput(ns("roi_poly_count")))
                ),

                actionButton(ns("roi_clear_selection"), "Clear current selection",
                             icon = icon("xmark"),
                             class = "btn-outline-secondary btn-sm w-100 mb-2"),

                textInput(ns("roi_name"), "ROI label", value = "roi"),
                actionButton(ns("roi_add"), "Add selection to ROI",
                             icon = icon("plus"), class = "btn-success w-100 mb-2"),
                actionButton(ns("roi_remove"), "Remove selection from ROI",
                             icon = icon("minus"), class = "btn-outline-danger w-100 mb-2"),
                actionButton(ns("roi_clear"), "Clear ROI",
                             icon = icon("eraser"), class = "btn-outline-secondary w-100"),
                div(class = "text-muted small mt-2", uiOutput(ns("roi_count")))
              )
            )
          ),

          column(8,
            bslib::card(
              bslib::card_header("Status"),
              bslib::card_body(uiOutput(ns("roi_status")))
            ),
            bslib::card(
              bslib::card_header("ROI map"),
              bslib::card_body(
                shinycssloaders::withSpinner(
                  plotOutput(ns("roi_plot"), height = "60vh",
                             click = clickOpts(id = ns("roi_click")),
                             brush = brushOpts(id = ns("roi_brush"), resetOnNew = FALSE))
                )
              )
            ),
            bslib::card(
              bslib::card_header("Export"),
              bslib::card_body(
                fluidRow(
                  column(6, downloadButton(ns("roi_export_csv"), "Download ROI (.csv)",
                                           class = "btn-outline-primary w-100")),
                  column(6, uiOutput(ns("save_roi_ui")))
                )
              )
            )
          )
        )
      )
    )
  )
}

spatialServer <- function(id, datasets, active_id) {
  moduleServer(id, function(input, output, session) {
    library(Cardinal)
    ns <- session$ns

    mse_data <- active_dataset(datasets, active_id)

    status <- reactiveVal(list(type = "info", message = "No segmentation run yet."))
    result <- reactiveVal(NULL)  # list(seg, table, plot_fn)

    output$active_dataset <- active_dataset_ui(datasets, active_id)

    observe({
      if (is.null(mse_data())) {
        shinyjs::disable("run_segmentation")
      } else {
        shinyjs::enable("run_segmentation")
      }
    })

    observeEvent(input$run_segmentation, {
      data <- mse_data()
      if (is.null(data)) {
        status(list(type = "danger", message = "No active dataset selected."))
        return()
      }

      # Cardinal's spatial engines fail on file-backed storage; materialize
      # into memory first (missing values imputed with 0).
      data <- tryCatch(
        with_progress_stages("Loading spectra into memory",
                     ensure_in_memory(data)),
        error = function(e) {
          status(list(type = "danger",
                      message = sprintf("Could not load spectra into memory: %s",
                                        conditionMessage(e))))
          NULL
        }
      )
      if (is.null(data)) return()

      seg <- with_progress_stages("Running spatial segmentation", {
        tryCatch(
          run_spatial_segmentation(
            data,
            r = if (is.null(input$seg_r)) 1 else input$seg_r,
            k = input$seg_k,
            s = if (is.null(input$seg_s)) 0 else input$seg_s,
            weights = if (is.null(input$seg_weights)) "gaussian" else input$seg_weights
          ),
          error = function(e) {
            status(list(type = "danger", message = conditionMessage(e)))
            NULL
          }
        )
      })
      if (is.null(seg)) return()

      result(list(
        kind = "segmentation",
        seg = seg,
        table = segmentation_top_features(seg, n = 50),
        plot_fn = function() plot_segmentation(seg)
      ))
      status(list(type = "success",
                  message = sprintf("Segmentation complete: %d segments over %d pixels (sizes: %s).",
                                    seg$n_segments, length(seg$segment),
                                    paste(as.integer(seg$sizes), collapse = ", "))))
    })

    output$result_plot <- renderPlot({
      res <- result()
      validate(need(!is.null(res), "Run a segmentation to see results."))
      res$plot_fn()
    })

    output$result_table <- DT::renderDT({
      res <- result()
      validate(need(!is.null(res), "Run a segmentation to see results."))
      tbl <- res$table
      num_cols <- vapply(tbl, is.numeric, logical(1))
      DT::datatable(tbl, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatSignif(columns = names(tbl)[num_cols], digits = 5)
    })

    output$save_segments_ui <- renderUI({
      res <- result()
      if (is.null(res)) {
        return(div(class = "text-muted small",
                   "Run a segmentation to enable saving."))
      }
      actionButton(ns("save_segments"), "Save segments as dataset",
                   icon = icon("floppy-disk"), class = "btn-outline-success w-100")
    })

    observeEvent(input$save_segments, {
      res <- result()
      req(!is.null(res))
      save_dataset_modal(ns, "Save segments as a new dataset",
                         "spatial_segments", "save_segments_confirm")
    })

    observeEvent(input$save_segments_confirm, {
      removeModal()
      res <- result()
      data <- mse_data()
      req(!is.null(res), !is.null(data))
      ds <- datasets()
      new_id <- next_dataset_id(names(ds))
      record <- tryCatch(
        make_segmentation_record(data, res$seg, new_id,
                                 filename = input$save_filename, ext = input$save_ext),
        error = function(e) {
          status(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(record)) return()
      ds[[new_id]] <- record
      datasets(ds)
      active_id(new_id)
      status(list(type = "success",
                  message = sprintf("Segments saved as '%s'.", record$name)))
    })

    output$export_csv <- downloadHandler(
      filename = function() {
        paste0("spatial_segments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        res <- result()
        if (is.null(res)) {
          utils::write.csv(data.frame(), file, row.names = FALSE)
        } else {
          utils::write.csv(segmentation_pixel_table(res$seg), file, row.names = FALSE)
        }
      }
    )

    output$export_plot <- export_plot_handler(
      result, function() input$plot_format, "spatial_segmentation",
      empty_msg = "No segmentation result to export.")

    output$status <- render_status(status)

    dgmm_status <- reactiveVal(list(type = "info",
                                    message = "No spatial-DGMM run yet."))
    dgmm_result <- reactiveVal(NULL)  # list(dgmm, table, plot_fn)

    observe({
      if (is.null(mse_data())) {
        shinyjs::disable("run_dgmm")
      } else {
        shinyjs::enable("run_dgmm")
      }
    })

    observeEvent(input$run_dgmm, {
      data <- mse_data()
      if (is.null(data)) {
        dgmm_status(list(type = "danger", message = "No active dataset selected."))
        return()
      }

      # spatialDGMM fails on file-backed storage (sparse gaps -> NA);
      # materialize into memory first (missing values imputed with 0).
      data <- tryCatch(
        with_progress_stages("Loading spectra into memory",
                     ensure_in_memory(data)),
        error = function(e) {
          dgmm_status(list(type = "danger",
                           message = sprintf("Could not load spectra into memory: %s",
                                             conditionMessage(e))))
          NULL
        }
      )
      if (is.null(data)) return()

      # Feature index from nearest m/z match, defaulting to feature 1.
      mz <- Cardinal::mz(data)
      feature <- 1L
      if (!is.null(input$dgmm_mz) && !is.na(input$dgmm_mz) && length(mz) > 0) {
        feature <- which.min(abs(mz - input$dgmm_mz))
      }

      dgmm <- with_progress_stages("Running spatial-DGMM", {
        tryCatch(
          run_spatial_dgmm(
            data,
            r = if (is.null(input$dgmm_r)) 1 else input$dgmm_r,
            k = input$dgmm_k,
            feature = feature,
            weights = if (is.null(input$dgmm_weights)) "gaussian" else input$dgmm_weights
          ),
          error = function(e) {
            dgmm_status(list(type = "danger", message = conditionMessage(e)))
            NULL
          }
        )
      })
      if (is.null(dgmm)) return()

      dgmm_result(list(
        kind = "dgmm",
        dgmm = dgmm,
        table = dgmm_pixel_table(dgmm),
        plot_fn = function() plot_dgmm(dgmm)
      ))
      mz_lab <- if (is.na(dgmm$mz)) sprintf("feature %d", dgmm$feature) else sprintf("m/z %.4g", dgmm$mz)
      dgmm_status(list(type = "success",
                       message = sprintf("Spatial-DGMM complete: %s segmented into %d components over %d pixels (sizes: %s).",
                                         mz_lab, dgmm$n_segments, length(dgmm$segment),
                                         paste(as.integer(dgmm$sizes), collapse = ", "))))
    })

    output$dgmm_plot <- renderPlot({
      res <- dgmm_result()
      validate(need(!is.null(res), "Run a spatial-DGMM to see results."))
      res$plot_fn()
    })

    output$dgmm_table <- DT::renderDT({
      res <- dgmm_result()
      validate(need(!is.null(res), "Run a spatial-DGMM to see results."))
      tbl <- res$table
      num_cols <- vapply(tbl, is.numeric, logical(1))
      DT::datatable(tbl, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatSignif(columns = names(tbl)[num_cols], digits = 5)
    })

    output$save_dgmm_ui <- renderUI({
      res <- dgmm_result()
      if (is.null(res)) {
        return(div(class = "text-muted small",
                   "Run a spatial-DGMM to enable saving."))
      }
      actionButton(ns("save_dgmm"), "Save segments as dataset",
                   icon = icon("floppy-disk"), class = "btn-outline-success w-100")
    })

    observeEvent(input$save_dgmm, {
      res <- dgmm_result()
      req(!is.null(res))
      save_dataset_modal(ns, "Save spatial-DGMM segments as a new dataset",
                         "dgmm_segments", "save_dgmm_confirm",
                         filename_id = "dgmm_save_filename",
                         ext_id = "dgmm_save_ext")
    })

    observeEvent(input$save_dgmm_confirm, {
      removeModal()
      res <- dgmm_result()
      data <- mse_data()
      req(!is.null(res), !is.null(data))
      ds <- datasets()
      new_id <- next_dataset_id(names(ds))
      record <- tryCatch(
        make_dgmm_record(data, res$dgmm, new_id,
                         filename = input$dgmm_save_filename, ext = input$dgmm_save_ext),
        error = function(e) {
          dgmm_status(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(record)) return()
      ds[[new_id]] <- record
      datasets(ds)
      active_id(new_id)
      dgmm_status(list(type = "success",
                       message = sprintf("Segments saved as '%s'.", record$name)))
    })

    output$dgmm_export_csv <- downloadHandler(
      filename = function() {
        paste0("dgmm_segments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        res <- dgmm_result()
        if (is.null(res)) {
          utils::write.csv(data.frame(), file, row.names = FALSE)
        } else {
          utils::write.csv(res$table, file, row.names = FALSE)
        }
      }
    )

    output$dgmm_export_plot <- export_plot_handler(
      dgmm_result, function() input$dgmm_plot_format, "dgmm_segmentation",
      empty_msg = "No spatial-DGMM result to export.")

    output$dgmm_status <- render_status(dgmm_status)

    # Per-pixel logical selection vector (length = ncol(active dataset)).
    roi_selected <- reactiveVal(NULL)
    # In-progress free-form polygon vertices (data.frame x/y) or NULL.
    roi_poly <- reactiveVal(NULL)
    roi_status <- reactiveVal(list(type = "info",
                                   message = "Drag a rectangle on the map to select pixels."))

    # Coordinates of the active dataset (x/y per pixel).
    roi_coord <- reactive({
      data <- mse_data()
      if (is.null(data)) return(NULL)
      as.data.frame(Cardinal::coord(data))
    })

    # Reset selection and outline whenever the active dataset changes.
    observeEvent(mse_data(), {
      data <- mse_data()
      if (is.null(data)) {
        roi_selected(NULL)
      } else {
        roi_selected(rep(FALSE, ncol(data)))
      }
      roi_poly(NULL)
    }, ignoreNULL = FALSE)

    # Drop the outline and any brush rectangle when the shape mode changes.
    observeEvent(input$roi_sel_mode, {
      roi_poly(NULL)
      session$resetBrush(ns("roi_brush"))
    })

    # Clear the uncommitted selection; the committed ROI stays.
    observeEvent(input$roi_clear_selection, {
      roi_poly(NULL)
      session$resetBrush(ns("roi_brush"))
      roi_status(list(type = "info", message = "Current selection cleared."))
    })

    # Append a clicked point to the polygon outline (polygon mode only).
    observeEvent(input$roi_click, {
      if (!identical(input$roi_sel_mode, "poly")) return()
      click <- input$roi_click
      if (is.null(click$x) || is.null(click$y)) return()
      poly <- roi_poly()
      pt <- data.frame(x = click$x, y = click$y)
      roi_poly(if (is.null(poly)) pt else rbind(poly, pt))
    })

    # Undo the last polygon vertex.
    observeEvent(input$roi_poly_undo, {
      poly <- roi_poly()
      if (is.null(poly) || nrow(poly) == 0) return()
      poly <- poly[-nrow(poly), , drop = FALSE]
      roi_poly(if (nrow(poly) == 0) NULL else poly)
    })

    # Clear the in-progress polygon outline (leaves the committed ROI intact).
    observeEvent(input$roi_poly_clear, {
      roi_poly(NULL)
    })

    # Polygon vertex count indicator.
    output$roi_poly_count <- renderUI({
      poly <- roi_poly()
      nv <- if (is.null(poly)) 0L else nrow(poly)
      msg <- sprintf("%d point%s in outline.", nv, if (nv == 1) "" else "s")
      if (nv > 0 && nv < 3) msg <- paste(msg, "Add at least 3 to enclose a region.")
      msg
    })

    # Pixels in the current selection; NULL + status message when empty.
    current_selection_hit <- function(coord) {
      if (identical(input$roi_sel_mode, "poly")) {
        poly <- roi_poly()
        if (is.null(poly) || nrow(poly) < 3) {
          roi_status(list(type = "warning",
                          message = "Click at least 3 points on the map to trace an outline first."))
          return(NULL)
        }
        hit <- pixels_in_polygon(coord, poly)
      } else {
        hit <- pixels_in_brush(coord, input$roi_brush)
      }
      if (!any(hit)) {
        roi_status(list(type = "warning",
                        message = "No pixels in the current selection. Make a selection on the map first."))
        return(NULL)
      }
      hit
    }

    # Union the selected pixels into the ROI.
    observeEvent(input$roi_add, {
      coord <- roi_coord()
      if (is.null(coord)) {
        roi_status(list(type = "danger", message = "No active dataset selected."))
        return()
      }
      hit <- current_selection_hit(coord)
      if (is.null(hit)) return()
      sel <- roi_selected()
      if (is.null(sel) || length(sel) != nrow(coord)) sel <- rep(FALSE, nrow(coord))
      sel <- sel | hit
      roi_selected(sel)
      if (identical(input$roi_sel_mode, "poly")) roi_poly(NULL)
      roi_status(list(type = "info",
                      message = sprintf("%d of %d pixels in ROI.", sum(sel), length(sel))))
    })

    # Remove the selected pixels from the ROI.
    observeEvent(input$roi_remove, {
      coord <- roi_coord()
      if (is.null(coord)) return()
      hit <- current_selection_hit(coord)
      if (is.null(hit)) return()
      sel <- roi_selected()
      if (is.null(sel) || length(sel) != nrow(coord)) sel <- rep(FALSE, nrow(coord))
      sel <- sel & !hit
      roi_selected(sel)
      if (identical(input$roi_sel_mode, "poly")) roi_poly(NULL)
      roi_status(list(type = "info",
                      message = sprintf("%d of %d pixels in ROI.", sum(sel), length(sel))))
    })

    # Clear the whole ROI (committed pixels, in-progress outline, and brush).
    observeEvent(input$roi_clear, {
      coord <- roi_coord()
      roi_selected(if (is.null(coord)) NULL else rep(FALSE, nrow(coord)))
      roi_poly(NULL)
      session$resetBrush(ns("roi_brush"))
      roi_status(list(type = "info", message = "ROI cleared."))
    })

    # Selected-pixel count.
    output$roi_count <- renderUI({
      sel <- roi_selected()
      n_sel <- if (is.null(sel)) 0L else sum(sel)
      total <- if (is.null(sel)) 0L else length(sel)
      sprintf("%d of %d pixels in ROI.", n_sel, total)
    })

    # ROI map; in polygon mode the in-progress outline is overlaid.
    output$roi_plot <- renderPlot({
      coord <- roi_coord()
      validate(need(!is.null(coord), "No active dataset. Import or select a dataset first."))
      roi_name <- if (isTruthy(input$roi_name)) input$roi_name else "roi"
      poly <- if (identical(input$roi_sel_mode, "poly")) roi_poly() else NULL
      plot_roi(coord, roi_selected(), roi_name = roi_name, poly = poly)
    })

    # ROI status output.
    output$roi_status <- render_status(roi_status)

    # Save-ROI button appears once at least one pixel is selected.
    output$save_roi_ui <- renderUI({
      sel <- roi_selected()
      if (is.null(sel) || sum(sel) == 0) {
        return(div(class = "text-muted small",
                   "Select pixels to enable saving."))
      }
      actionButton(ns("save_roi"), "Save ROI as dataset",
                   icon = icon("floppy-disk"), class = "btn-outline-success w-100")
    })

    observeEvent(input$save_roi, {
      sel <- roi_selected()
      req(!is.null(sel), sum(sel) > 0)
      save_dataset_modal(ns, "Save ROI as a new dataset",
                         "roi_selection", "save_roi_confirm",
                         filename_id = "roi_save_filename",
                         ext_id = "roi_save_ext")
    })

    observeEvent(input$save_roi_confirm, {
      removeModal()
      sel <- roi_selected()
      data <- mse_data()
      req(!is.null(sel), !is.null(data), sum(sel) > 0)
      roi_name <- if (isTruthy(input$roi_name)) input$roi_name else "roi"
      ds <- datasets()
      new_id <- next_dataset_id(names(ds))
      record <- tryCatch(
        make_roi_record(data, sel, new_id, roi_name = roi_name, column = roi_name,
                        filename = input$roi_save_filename, ext = input$roi_save_ext),
        error = function(e) {
          roi_status(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(record)) return()
      ds[[new_id]] <- record
      datasets(ds)
      active_id(new_id)
      roi_status(list(type = "success",
                      message = sprintf("ROI saved as '%s'.", record$name)))
    })

    output$roi_export_csv <- downloadHandler(
      filename = function() {
        paste0("roi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        coord <- roi_coord()
        sel <- roi_selected()
        if (is.null(coord)) {
          utils::write.csv(data.frame(), file, row.names = FALSE)
        } else {
          roi_name <- if (isTruthy(input$roi_name)) input$roi_name else "roi"
          utils::write.csv(roi_pixel_table(coord, sel, roi_name), file, row.names = FALSE)
        }
      }
    )
  })
}
