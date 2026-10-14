# Cardinal NEST - preprocessing module: ordered pipeline builder and runner

library(shiny)

# Each entry maps to a Cardinal function and drives the dynamic UI; adding a
# step only requires editing this list. Param `type` is
# "select"|"numeric"|"text"; `show_if` shows a param conditionally.
preproc_step_defs <- list(
  normalize = list(
    label = "Normalization",
    fn = "normalize",
    group = "core",
    desc = "Scale each spectrum (total ion current, RMS, or to a reference).",
    params = list(
      list(id = "method", label = "Method", type = "select",
           choices = c("TIC" = "tic", "RMS" = "rms", "Reference" = "reference"),
           default = "tic"),
      list(id = "scale", label = "Scale (optional)", type = "numeric",
           default = NA, step = "any"),
      list(id = "ref", label = "Reference m/z value(s)", type = "text",
           default = "", placeholder = "e.g. 1234.5, 1567.8",
           show_if = list(input = "method", value = "reference"))
    )
  ),
  reduceBaseline = list(
    label = "Baseline Correction",
    fn = "reduceBaseline",
    group = "core",
    desc = "Estimate and subtract the baseline from each spectrum.",
    params = list(
      list(id = "method", label = "Method", type = "select",
           choices = c("Local minima" = "locmin", "Convex hull" = "hull",
                       "SNIP" = "snip", "Median" = "median"),
           default = "locmin")
    )
  ),
  peakPick = list(
    label = "Peak Picking",
    fn = "peakPick",
    group = "core",
    desc = "Detect peaks in each spectrum using a noise-based threshold.",
    params = list(
      list(id = "method", label = "Method", type = "select",
           choices = c("Derivative" = "diff", "SD" = "sd", "MAD" = "mad",
                       "Quantile" = "quantile", "Filter" = "filter", "CWT" = "cwt"),
           default = "diff"),
      list(id = "SNR", label = "Signal-to-noise ratio", type = "numeric",
           default = 2, step = "any"),
      list(id = "type", label = "Peak value", type = "select",
           choices = c("Height" = "height", "Area" = "area"),
           default = "height"),
      list(id = "tolerance", label = "Tolerance (optional)", type = "numeric",
           default = NA, step = "any"),
      list(id = "units", label = "Tolerance Units", type = "select",
           choices = c("ppm" = "ppm", "m/z" = "mz"), default = "ppm")
    )
  ),
  peakAlign = list(
    label = "Peak Alignment",
    fn = "peakAlign",
    group = "core",
    desc = "Align detected peaks across spectra to a common set of m/z values.",
    params = list(
      list(id = "tolerance", label = "Tolerance (optional)", type = "numeric",
           default = NA, step = "any"),
      list(id = "units", label = "Tolerance Units", type = "select",
           choices = c("ppm" = "ppm", "m/z" = "mz"), default = "ppm")
    )
  ),
  smooth = list(
    label = "Smoothing",
    fn = "smooth",
    group = "optional",
    desc = "Reduce noise by smoothing each spectrum.",
    params = list(
      list(id = "method", label = "Method", type = "select",
           choices = c("Gaussian" = "gaussian", "Bilateral" = "bilateral",
                       "Adaptive" = "adaptive", "Guided" = "guide",
                       "Peak-aware guided" = "pag",
                       "Savitzky-Golay" = "sgolay", "Moving average" = "ma"),
           default = "gaussian"),
      list(id = "width", label = "Window width", type = "numeric",
           default = 11, step = 1)
    )
  ),
  recalibrate = list(
    label = "Recalibration",
    fn = "recalibrate",
    group = "optional",
    desc = "Correct systematic m/z drift against reference peaks.",
    params = list(
      list(id = "method", label = "Method", type = "select",
           choices = c("Local maxima" = "locmax", "DTW" = "dtw", "COW" = "cow"),
           default = "locmax"),
      list(id = "tolerance", label = "Tolerance (optional)", type = "numeric",
           default = NA, step = "any"),
      list(id = "units", label = "Tolerance Units", type = "select",
           choices = c("ppm" = "ppm", "m/z" = "mz"), default = "ppm")
    )
  ),
  peakFilter = list(
    label = "Peak Filtering",
    fn = "peakFilter",
    group = "optional",
    desc = paste("Remove low-frequency peaks by minimum occurrence.",
                 "Note: peakFilter is deprecated in this Cardinal version;",
                 "filtering will be applied via subsetFeatures on peak frequency."),
    params = list(
      list(id = "freq_min", label = "Minimum frequency (0-1)", type = "numeric",
           default = 0.01, step = "any")
    )
  ),
  bin = list(
    label = "Binning",
    fn = "bin",
    group = "optional",
    desc = "Bin spectra to a common m/z axis at a chosen resolution.",
    params = list(
      list(id = "resolution", label = "Resolution", type = "numeric",
           default = 1, step = "any"),
      list(id = "units", label = "Resolution Units", type = "select",
           choices = c("ppm" = "ppm", "m/z" = "mz"), default = "mz"),
      list(id = "method", label = "Summary", type = "select",
           choices = c("Sum" = "sum", "Mean" = "mean",
                       "Max" = "max", "Min" = "min"),
           default = "sum")
    )
  )
)

# Choices for the "add step" selector.
preproc_step_choices <- function() {
  steps <- names(preproc_step_defs)
  labels <- vapply(preproc_step_defs, function(d) d$label, character(1))
  setNames(steps, labels)
}

preprocessingUI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
             h3("Preprocessing"),
             p("Build a custom pipeline of Cardinal preprocessing steps and run it on the active dataset."),
             uiOutput(ns("active_dataset"))
      )
    ),

    fluidRow(
      column(4,
        bslib::card(
          bslib::card_header("Add a step"),
          bslib::card_body(
            selectInput(ns("step_type"), "Preprocessing step",
                        choices = preproc_step_choices()),
            div(class = "text-muted small mb-2", uiOutput(ns("step_type_desc"))),
            actionButton(ns("add_step"), "Add to pipeline",
                         icon = icon("plus"), class = "btn-primary w-100")
          )
        ),
        bslib::card(
          bslib::card_header("Run"),
          bslib::card_body(
            actionButton(ns("run_pipeline"), "Run pipeline",
                         icon = icon("play"), class = "btn-success w-100 mb-2"),
            actionButton(ns("clear_pipeline"), "Clear pipeline",
                         icon = icon("trash"), class = "btn-outline-danger w-100"),
            checkboxInput(ns("save_as_new"),
                          "Save result as a new dataset", value = TRUE)
          )
        )
      ),

      column(8,
        bslib::card(
          bslib::card_header("Pipeline"),
          bslib::card_body(
            uiOutput(ns("pipeline_steps"))
          )
        ),
        bslib::card(
          bslib::card_header("Status"),
          bslib::card_body(
            uiOutput(ns("run_status"))
          )
        )
      )
    )
  )
}

# One parameter control for a step instance (uid) and param def.
render_step_param <- function(ns, uid, param) {
  input_id <- ns(paste0("p_", uid, "_", param$id))
  control <- switch(param$type,
    "select" = selectInput(input_id, param$label,
                           choices = param$choices, selected = param$default),
    "numeric" = numericInput(input_id, param$label,
                             value = param$default,
                             step = if (!is.null(param$step)) param$step else "any"),
    "text" = textInput(input_id, param$label,
                       value = if (!is.null(param$default)) as.character(param$default) else "",
                       placeholder = if (!is.null(param$placeholder)) param$placeholder else ""),
    NULL
  )

  if (is.null(control)) return(NULL)

  if (!is.null(param$show_if)) {
    trigger_id <- ns(paste0("p_", uid, "_", param$show_if$input))
    condition <- sprintf("input['%s'] == '%s'", trigger_id, param$show_if$value)
    conditionalPanel(condition = condition, control)
  } else {
    control
  }
}

# Step card with reorder/remove/collapse controls.
render_step_card <- function(ns, step, index, n_steps, collapsed = FALSE) {
  def <- preproc_step_defs[[step$type]]
  uid <- step$uid

  # Action + uid encoded into a single namespaced input event.
  action_btn <- function(action, icon_name, title, disabled = FALSE) {
    tags$button(
      type = "button",
      class = paste("btn btn-sm btn-outline-secondary", if (disabled) "disabled" else ""),
      title = title,
      onclick = sprintf(
        "Shiny.setInputValue('%s', {uid: '%s', nonce: Math.random()}, {priority: 'event'})",
        ns(action), uid
      ),
      icon(icon_name)
    )
  }

  bslib::card(
    class = "mb-2",
    bslib::card_header(
      div(class = "d-flex justify-content-between align-items-center",
        span(strong(sprintf("%d. %s", index, def$label))),
        div(class = "btn-group",
          action_btn("collapse_step",
                     if (collapsed) "chevron-right" else "chevron-down",
                     if (collapsed) "Expand" else "Collapse"),
          action_btn("move_up", "arrow-up", "Move up", disabled = index == 1),
          action_btn("move_down", "arrow-down", "Move down", disabled = index == n_steps),
          action_btn("remove_step", "times", "Remove")
        )
      )
    ),
    tags$div(id = ns(paste0("step_body_", uid)),
             style = if (collapsed) "display: none;" else "",
      bslib::card_body(
        p(class = "text-muted small", def$desc),
        fluidRow(
          lapply(def$params, function(param) {
            column(4, render_step_param(ns, uid, param))
          })
        )
      )
    )
  )
}

# Run `expr` while capturing warnings (re-emitted) tagged with `label`.
# Returns list(value, warnings).
capture_pipeline_warnings <- function(expr, label = "") {
  captured <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      wm <- conditionMessage(w)
      captured[[length(captured) + 1L]] <<-
        if (nzchar(label)) sprintf("%s: %s", label, wm) else wm
    }
  )
  list(value = value, warnings = captured)
}

run_preproc_pipeline <- function(data, specs, progress = NULL) {
  if (is.null(data)) stop("No dataset supplied to the pipeline.", call. = FALSE)

  warnings <- character()
  result <- data

  # `progress` (optional): function(frac, detail); bar spans n steps + the
  # final process() call.
  n_steps <- length(specs)
  report <- function(done, detail) {
    if (!is.null(progress)) progress(done / (n_steps + 1L), detail)
  }

  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    def <- preproc_step_defs[[spec$type]]
    if (is.null(def)) {
      stop(sprintf("Step %d: unknown preprocessing step '%s'.", i, spec$type),
           call. = FALSE)
    }
    label <- def$label
    report(i - 1L, sprintf("Step %d/%d: %s", i, n_steps, label))

    # peakFilter is deprecated in this Cardinal version; skip it.
    if (spec$type == "peakFilter") {
      warnings <- c(warnings, "peakFilter is deprecated and was skipped.")
      next
    }

    # Peak picking already-centroided data is meaningless and can overflow
    # on large datasets; skip with a warning.
    if (spec$type == "peakPick" && isTRUE(tryCatch(centroided(result),
                                                   error = function(e) FALSE))) {
      warnings <- c(warnings,
                    "Peak Picking skipped: dataset is already centroided.")
      next
    }

    fn <- tryCatch(get(def$fn, envir = asNamespace("Cardinal")),
                   error = function(e) NULL)
    if (is.null(fn)) {
      stop(sprintf("Step %d (%s): Cardinal function '%s' not found.",
                   i, label, def$fn), call. = FALSE)
    }

    step_label <- sprintf("Step %d (%s)", i, label)
    run <- tryCatch(
      capture_pipeline_warnings(do.call(fn, c(list(result), spec$args)), step_label),
      error = function(e) {
        stop(sprintf("Step %d (%s): %s", i, label, conditionMessage(e)),
             call. = FALSE)
      }
    )
    result <- run$value
    warnings <- c(warnings, run$warnings)
  }

  # Flush queued transformations (normalize/smooth/reduceBaseline/peakPick).
  report(n_steps, "Executing queued steps")
  run <- tryCatch(
    capture_pipeline_warnings(Cardinal::process(result), "Executing queued steps"),
    error = function(e) {
      stop(sprintf("Executing queued steps failed: %s", conditionMessage(e)),
           call. = FALSE)
    }
  )
  result <- run$value
  warnings <- c(warnings, run$warnings)
  report(n_steps + 1L, "Complete")

  list(data = result, warnings = warnings)
}

# Generate a dataset id that does not collide with any existing id.
next_dataset_id <- function(existing_ids) {
  i <- length(existing_ids) + 1
  repeat {
    candidate <- paste0("ds_", i)
    if (!candidate %in% existing_ids) return(candidate)
    i <- i + 1
  }
}

# Dataset record for a preprocessed result, matching the Import tab's
# structure; `ext` is appended to the display name.
make_preproc_record <- function(result, id, step_labels, filename = NULL, ext = NULL) {
  nm <- if (!is.null(filename)) trimws(filename) else ""
  if (!nzchar(nm)) nm <- "preprocessed"
  if (!is.null(ext) && nzchar(ext)) nm <- paste0(nm, ".", ext)

  # `centroided` comes from ProtGenerics via Cardinal's search path; it is
  # not re-exported so it must stay unqualified.
  summary <- list(
    n_pixels = ncol(result),
    n_features = nrow(result),
    n_runs = length(Cardinal::runNames(result)),
    mz_range = range(Cardinal::mz(result), na.rm = TRUE),
    centroided = centroided(result),
    class = class(result)[1]
  )

  list(
    id = id,
    name = nm,
    format = "preprocessed",
    extensions = paste(step_labels, collapse = " -> "),
    files = "",
    imported_at = Sys.time(),
    summary = summary,
    data = result
  )
}

preprocessingServer <- function(id, datasets, active_id, perf_log = NULL) {
  moduleServer(id, function(input, output, session) {
    library(Cardinal)
    ns <- session$ns

    mse_data <- active_dataset(datasets, active_id)

    # Ordered pipeline: list of steps, each list(uid, type).
    pipeline <- reactiveVal(list())
    step_counter <- reactiveVal(0)
    run_status <- reactiveVal(list(type = "info", message = "No pipeline run yet."))
    collapsed_uids <- reactiveVal(character())

    output$active_dataset <- active_dataset_ui(datasets, active_id)

    output$step_type_desc <- renderUI({
      req(input$step_type)
      def <- preproc_step_defs[[input$step_type]]
      if (is.null(def)) return(NULL)
      def$desc
    })

    observeEvent(input$add_step, {
      req(input$step_type)
      new_count <- step_counter() + 1
      step_counter(new_count)
      uid <- paste0("s", new_count)
      current <- pipeline()
      current[[length(current) + 1]] <- list(uid = uid, type = input$step_type)
      pipeline(current)
    })

    observeEvent(input$clear_pipeline, {
      pipeline(list())
    })

    observeEvent(input$remove_step, {
      uid <- input$remove_step$uid
      current <- pipeline()
      keep <- vapply(current, function(s) s$uid != uid, logical(1))
      pipeline(current[keep])
    })

    observeEvent(input$move_up, {
      uid <- input$move_up$uid
      current <- pipeline()
      idx <- which(vapply(current, function(s) s$uid == uid, logical(1)))
      if (length(idx) == 1 && idx > 1) {
        current[c(idx - 1, idx)] <- current[c(idx, idx - 1)]
        pipeline(current)
      }
    })

    observeEvent(input$move_down, {
      uid <- input$move_down$uid
      current <- pipeline()
      idx <- which(vapply(current, function(s) s$uid == uid, logical(1)))
      if (length(idx) == 1 && idx < length(current)) {
        current[c(idx, idx + 1)] <- current[c(idx + 1, idx)]
        pipeline(current)
      }
    })

    observeEvent(input$collapse_step, {
      uid <- input$collapse_step$uid
      current <- collapsed_uids()
      if (uid %in% current) {
        collapsed_uids(setdiff(current, uid))
      } else {
        collapsed_uids(union(current, uid))
      }
    })

    output$pipeline_steps <- renderUI({
      steps <- pipeline()
      if (length(steps) == 0) {
        return(div(class = "text-muted",
                   "Your pipeline is empty. Add a step from the panel on the left."))
      }
      n <- length(steps)
      collapsed <- collapsed_uids()
      tagList(
        lapply(seq_along(steps), function(i) {
          step <- steps[[i]]
          is_collapsed <- step$uid %in% collapsed
          render_step_card(ns, step, i, n, collapsed = is_collapsed)
        })
      )
    })

    observe({
      if (is.null(mse_data()) || length(pipeline()) == 0) {
        shinyjs::disable("run_pipeline")
      } else {
        shinyjs::enable("run_pipeline")
      }
    })

    # Collect arguments for a single step from its UI inputs.
    build_step_args <- function(step) {
      def <- preproc_step_defs[[step$type]]
      args <- list()
      for (param in def$params) {
        input_id <- paste0("p_", step$uid, "_", param$id)
        val <- input[[input_id]]
        if (is.null(val)) next
        if (param$type == "numeric") {
          if (length(val) == 1 && is.na(val)) next
          args[[param$id]] <- as.numeric(val)
        } else if (param$type == "text") {
          if (!is.character(val) || !nzchar(as.character(val))) next
          args[[param$id]] <- as.character(val)
        } else {
          args[[param$id]] <- val
        }
      }

      # normalize ref: only used for method = "reference"
      if (step$type == "normalize") {
        if (is.null(args$method) || args$method != "reference") {
          args$ref <- NULL
        } else if (!is.null(args$ref)) {
          nums <- suppressWarnings(as.numeric(unlist(strsplit(args$ref, "[,;[:space:]]+"))))
          nums <- nums[!is.na(nums)]
          args$ref <- if (length(nums) == 0) NULL else nums
        }
        if (!is.null(args$scale) && length(args$scale) == 1 && is.na(args$scale)) {
          args$scale <- NULL
        }
      }

      # units are only meaningful when a numeric partner is supplied
      if (step$type %in% c("peakPick", "peakAlign", "recalibrate")) {
        if (is.null(args$tolerance) || (length(args$tolerance) == 1 && is.na(args$tolerance))) {
          args$tolerance <- NULL
          args$units <- NULL
        }
      }

      args
    }

    # Execute the pipeline and optionally save the result as a new dataset.
    execute_pipeline <- function(save = FALSE, filename = NULL, ext = NULL) {
      data <- mse_data()
      if (is.null(data)) {
        run_status(list(type = "danger", message = "No active dataset selected."))
        return()
      }
      steps <- pipeline()
      if (length(steps) == 0) {
        run_status(list(type = "warning", message = "No pipeline steps to run."))
        return()
      }

      # Shiny-independent step specs from the configured UI inputs.
      specs <- lapply(steps, function(step) {
        list(type = step$type, args = build_step_args(step))
      })

      pipe_t0 <- proc.time()[["elapsed"]]
      run <- tryCatch(
        withProgress(message = "Running preprocessing pipeline", value = 0, {
          run_preproc_pipeline(
            data, specs,
            progress = function(frac, detail) {
              setProgress(value = frac, detail = detail)
            })
        }),
        error = function(e) {
          run_status(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(run)) return()
      if (!is.null(perf_log)) {
        perf_log(perf_log_append(perf_log(), "preprocess",
                                 sprintf("pipeline (%d step(s))", length(steps)),
                                 proc.time()[["elapsed"]] - pipe_t0))
      }

      result <- run$data
      warnings <- run$warnings

      saved_name <- NULL
      if (save) {
        ds <- datasets()
        new_id <- next_dataset_id(names(ds))
        step_labels <- vapply(steps, function(s) preproc_step_defs[[s$type]]$label,
                              character(1))
        record <- make_preproc_record(result, new_id, step_labels,
                                       filename = filename, ext = ext)
        ds[[new_id]] <- record
        datasets(ds)
        active_id(new_id)
        saved_name <- record$name
      }

      msg <- paste0("Pipeline completed (", length(steps), " step(s))")
      if (!is.null(saved_name)) {
        msg <- paste0(msg, ". Saved as '", saved_name, "'.")
      } else {
        msg <- paste0(msg, ".")
      }
      # Warnings do not mean failure; they render as a callout below.
      run_status(list(type = "success", message = msg, warnings = warnings))
    }

    observeEvent(input$run_pipeline, {
      if (isTRUE(input$save_as_new)) {
        save_dataset_modal(ns, "Save preprocessed dataset",
                           "preprocessed_dataset", "save_run",
                           confirm_label = "Run and save")
      } else {
        execute_pipeline(save = FALSE)
      }
    })

    observeEvent(input$save_run, {
      removeModal()
      execute_pipeline(save = TRUE,
                       filename = input$save_filename,
                       ext = input$save_ext)
    })

    output$run_status <- renderUI({
      st <- run_status()
      cls <- status_alert_class(st$type)
      warns <- st$warnings
      # Warnings render as an amber callout regardless of status colour.
      warning_ui <- if (!is.null(warns) && length(warns) > 0) {
        div(class = "alert alert-warning mt-2 mb-0",
            div(strong(sprintf("%d pipeline warning(s):", length(warns)))),
            tags$ul(class = "mb-0",
                    lapply(warns, function(w) tags$li(w))))
      } else {
        NULL
      }
      tagList(
        div(class = paste("alert", cls, if (!is.null(warning_ui)) "mb-2" else ""),
            st$message),
        warning_ui
      )
    })

    # Current pipeline as Shiny-independent specs, consumed by the Workflow
    # tab for save/export/batch.
    current_specs <- reactive({
      steps <- pipeline()
      lapply(steps, function(step) {
        list(type = step$type, args = build_step_args(step))
      })
    })

    return(list(specs = current_specs))
  })
}
