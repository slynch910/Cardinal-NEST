# Cardinal NEST - performance module
# BiocParallel backend, dataset memory and conversion, result cache,
# global seed, and a timing log.

# MulticoreParam forks, unavailable on Windows.
available_backends <- function(os = .Platform$OS.type) {
  backends <- c("Serial (single core)" = "serial",
                "Snow (multiple processes)" = "snow")
  if (!identical(os, "windows")) {
    backends <- c(backends, "Multicore (fork)" = "multicore")
  }
  backends
}

# `seed` becomes RNGseed so bplapply streams match across backends.
make_bpparam <- function(backend = "serial", workers = 1L, seed = NULL) {
  backend <- match.arg(backend, c("serial", "snow", "multicore"))
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) workers <- 1L
  seed <- if (is.null(seed)) NA_integer_ else parse_seed(seed)
  if (is.na(seed)) seed <- NULL
  switch(backend,
         serial = BiocParallel::SerialParam(progressbar = FALSE, RNGseed = seed),
         snow = BiocParallel::SnowParam(workers = workers, progressbar = FALSE,
                                        RNGseed = seed),
         multicore = BiocParallel::MulticoreParam(workers = workers,
                                                  progressbar = FALSE,
                                                  RNGseed = seed))
}

describe_bpparam <- function(param) {
  cls <- class(param)[1]
  workers <- tryCatch(BiocParallel::bpnworkers(param),
                      error = function(e) NA_integer_)
  if (is.na(workers)) return(cls)
  seed <- tryCatch(BiocParallel::bpRNGseed(param), error = function(e) NULL)
  sprintf("%s, %d worker%s%s", cls, workers, if (workers == 1L) "" else "s",
          if (is.null(seed)) "" else sprintf(", RNG seed %s", seed))
}

# Distinct PIDs (none equal to Sys.getpid()) confirm real workers.
backend_worker_pids <- function(param, tasks = NULL) {
  if (is.null(tasks)) tasks <- max(2L, 2L * BiocParallel::bpnworkers(param))
  unlist(BiocParallel::bplapply(seq_len(tasks), function(i) Sys.getpid(),
                                BPPARAM = param))
}

# Worker PIDs plus serial-vs-`param` summarizeFeatures() timing when `data`
# is given; workers start beforehand so startup cost is excluded.
test_backend <- function(param, data = NULL) {
  if (!BiocParallel::bpisup(param)) {
    BiocParallel::bpstart(param)
    on.exit(BiocParallel::bpstop(param), add = TRUE)
  }
  pids <- unique(backend_worker_pids(param))
  serial_s <- parallel_s <- NA_real_
  if (!is.null(data)) {
    summarize <- function(p) {
      with_timing(Cardinal::summarizeFeatures(data, stat = "mean",
                                              verbose = FALSE, BPPARAM = p))$seconds
    }
    serial_s <- summarize(BiocParallel::SerialParam(progressbar = FALSE))
    parallel_s <- summarize(param)
  }
  list(backend = describe_bpparam(param),
       parallel = !methods::is(param, "SerialParam"),
       main_pid = Sys.getpid(),
       pids = pids,
       n_workers_used = length(pids),
       serial_seconds = serial_s,
       parallel_seconds = parallel_s,
       speedup = if (is.na(parallel_s)) NA_real_ else serial_s / max(parallel_s, 1e-6))
}

backend_test_summary <- function(res) {
  workers_ok <- !res$parallel ||
    (res$n_workers_used > 1L && !(res$main_pid %in% res$pids))
  msg <- sprintf("%s: %d worker process%s used (PID %s).", res$backend,
                 res$n_workers_used, if (res$n_workers_used == 1L) "" else "es",
                 paste(res$pids, collapse = ", "))
  if (!workers_ok) {
    msg <- paste(msg, "Expected more than one worker process separate from the app.")
  }
  if (is.na(res$speedup)) {
    msg <- paste(msg, "Load a dataset to also time summarizeFeatures.")
  } else {
    msg <- paste(msg, sprintf("summarizeFeatures: serial %.2fs, this backend %.2fs (%.2fx speedup).",
                              res$serial_seconds, res$parallel_seconds, res$speedup))
  }
  type <- if (!workers_ok) "warning" else if (is.na(res$speedup)) "info" else "success"
  list(type = type, message = msg)
}

# Per-dataset memory table: name, dims, estimated size, storage mode.
dataset_memory_table <- function(datasets, size_fn = estimate_spectra_bytes) {
  if (length(datasets) == 0) {
    return(data.frame(Dataset = character(0), Pixels = integer(0),
                      Features = integer(0), `Est. size` = character(0),
                      Storage = character(0), check.names = FALSE))
  }
  rows <- lapply(datasets, function(rec) {
    d <- rec$data
    bytes <- tryCatch(size_fn(d), error = function(e) NA_real_)
    data.frame(
      Dataset = rec$name %||% rec$id,
      Pixels = tryCatch(ncol(d), error = function(e) NA_integer_),
      Features = tryCatch(nrow(d), error = function(e) NA_integer_),
      `Est. size` = format_bytes(bytes),
      Storage = storage_label(d),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

performanceUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
             h3("Performance Tuning"),
             p(paste("Configure parallel processing, manage dataset memory, and",
                     "monitor operation timings."))
      )
    ),
    fluidRow(
      column(5,
             bslib::card(
               bslib::card_header("Parallel processing"),
               bslib::card_body(
                 selectInput(ns("backend"), "Backend",
                             choices = available_backends()),
                 conditionalPanel(
                   condition = sprintf("input['%s'] != 'serial'", ns("backend")),
                   numericInput(ns("workers"), "Workers", value = 2,
                                min = 1, max = 64, step = 1)
                 ),
                 div(class = "d-flex gap-2",
                     actionButton(ns("apply_backend"), "Apply backend",
                                  icon = icon("microchip"), class = "btn-primary"),
                     actionButton(ns("test_backend"), "Test backend",
                                  icon = icon("stethoscope"),
                                  class = "btn-outline-secondary")),
                 div(class = "text-muted small mt-2",
                     paste("The global random seed below is applied as the backend's",
                           "RNG seed when you click Apply backend.")),
                 uiOutput(ns("backend_status")),
                 uiOutput(ns("backend_test"))
               )
             ),
             bslib::card(
               bslib::card_header("Result cache"),
               bslib::card_body(
                 checkboxInput(ns("cache_enabled"),
                               "Cache analysis results (re-run identical inputs instantly)",
                               value = TRUE),
                 div(class = "text-muted small mb-2", textOutput(ns("cache_size"))),
                 actionButton(ns("clear_cache"), "Clear cache",
                              icon = icon("trash"), class = "btn-outline-secondary btn-sm")
               )
             ),
             bslib::card(
               bslib::card_header("Reproducibility"),
               bslib::card_body(
                 seed_input(ns, "global_seed", label = "Global random seed"),
                 div(class = "text-muted small",
                     paste("Used by PCA, NMF, k-means, cross-validation, and spatial",
                           "segmentation whenever their own seed is left blank."))
               )
             )
      ),
      column(7,
             bslib::card(
               bslib::card_header("Dataset memory"),
               bslib::card_body(
                 tableOutput(ns("memory_table")),
                 uiOutput(ns("memory_budget")),
                 uiOutput(ns("convert_ui")),
                 uiOutput(ns("convert_status"))
               )
             ),
             bslib::card(
               bslib::card_header("Operation timings"),
               bslib::card_body(
                 p(class = "text-muted small",
                   "Elapsed wall-clock time of recent imports, preprocessing, and analyses."),
                 tableOutput(ns("timing_table"))
               )
             )
      )
    )
  )
}

performanceServer <- function(id, datasets, active_id,
                              perf_log = reactiveVal(NULL),
                              perf_state = reactiveValues(cache_enabled = TRUE),
                              analysis_cache = new_result_cache()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Init-time read outside a reactive consumer, hence isolate().
    updateNumericInput(session, "global_seed",
                       value = shiny::isolate(session_global_seed(perf_state)))
    observe_seed_generators(input, session, "global_seed")
    observeEvent(input$global_seed, {
      s <- parse_seed(input$global_seed)
      if (!is.na(s)) perf_state$global_seed <- s
    })

    backend_msg <- reactiveVal(NULL)
    backend_test_msg <- reactiveVal(NULL)
    convert_msg <- reactiveVal(NULL)

    observeEvent(input$cache_enabled, {
      perf_state$cache_enabled <- isTRUE(input$cache_enabled)
    }, ignoreInit = FALSE)

    observeEvent(input$apply_backend, {
      param <- tryCatch(
        make_bpparam(input$backend, input$workers,
                     seed = session_global_seed(perf_state)),
        error = function(e) {
          backend_msg(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(param)) return()
      ok <- tryCatch({
        Cardinal::setCardinalBPPARAM(param)
        TRUE
      }, error = function(e) {
        backend_msg(list(type = "danger",
                         message = paste("Failed to set backend:", conditionMessage(e))))
        FALSE
      })
      if (ok) {
        backend_msg(list(type = "success",
                         message = paste("Backend set to", describe_bpparam(param))))
      }
    })

    output$backend_status <- renderUI({
      msg <- backend_msg()
      current <- tryCatch(describe_bpparam(Cardinal::getCardinalBPPARAM()),
                          error = function(e) "unknown")
      tagList(
        div(class = "text-muted small mt-2",
            strong("Current: "), current),
        if (!is.null(msg)) status_alert(msg, extra_class = "mt-2")
      )
    })

    observeEvent(input$test_backend, {
      ds <- datasets()
      aid <- active_id()
      data <- if (!is.null(aid) && !is.null(ds[[aid]])) ds[[aid]]$data else NULL
      timed <- tryCatch(
        withProgress(message = "Testing parallel backend", value = 0.2,
                     with_timing(test_backend(Cardinal::getCardinalBPPARAM(), data))),
        error = function(e) {
          backend_test_msg(list(type = "danger",
                                message = paste("Backend test failed:", conditionMessage(e))))
          NULL
        }
      )
      if (is.null(timed)) return()
      perf_log(perf_log_append(perf_log(), "performance", "test backend",
                               timed$seconds, note = timed$value$backend))
      backend_test_msg(backend_test_summary(timed$value))
    })

    output$backend_test <- renderUI({
      msg <- backend_test_msg()
      if (is.null(msg)) return(NULL)
      status_alert(msg, extra_class = "mt-2")
    })

    output$cache_size <- renderText({
      # Depend on perf_state$cache_version so this re-renders when any module
      # mutates the shared cache (Analyze stores, clear_cache empties).
      perf_state$cache_version
      sprintf("Cached results: %d", analysis_cache$size())
    })

    observeEvent(input$clear_cache, {
      analysis_cache$clear()
      perf_state$cache_version <-
        (if (is.null(perf_state$cache_version)) 0L else perf_state$cache_version) + 1L
    })

    output$memory_table <- renderTable({
      dataset_memory_table(datasets())
    })

    output$memory_budget <- renderUI({
      datasets()
      avail <- available_memory_bytes()
      div(class = "text-muted small mb-2",
          strong("Available memory: "),
          if (is.na(avail)) "unknown" else format_bytes(avail),
          sprintf(" (conversions warn above %.0f%% and are blocked above %.0f%% of it)",
                  100 * getOption("cardinalnest.memory_warn_frac", 0.5),
                  100 * getOption("cardinalnest.memory_block_frac", 0.8)))
    })

    output$convert_ui <- renderUI({
      ds <- datasets()
      aid <- active_id()
      if (is.null(aid) || is.null(ds[[aid]])) return(NULL)
      if (identical(spectra_storage(ds[[aid]]$data), "memory")) {
        return(div(class = "text-muted small",
                   "Active dataset is already in memory."))
      }
      actionButton(ns("convert_active"),
                   "Convert active dataset to in-memory",
                   icon = icon("bolt"), class = "btn-outline-primary btn-sm")
    })

    observeEvent(input$convert_active, {
      ds <- datasets()
      aid <- active_id()
      if (is.null(aid) || is.null(ds[[aid]])) return()
      pf <- memory_preflight(
        estimate_spectra_bytes(ds[[aid]]$data),
        action = sprintf("Converting '%s' to in-memory", ds[[aid]]$name),
        alternative = "Keep it file-backed; preprocessing and analyses run on file-backed data.")
      if (identical(pf$level, "block")) {
        convert_msg(list(type = "danger", message = pf$message))
        return()
      }
      if (identical(pf$level, "warn")) {
        showModal(modalDialog(
          title = "Low memory",
          pf$message,
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("confirm_convert"), "Convert anyway",
                         class = "btn-warning")
          )
        ))
        return()
      }
      convert_active_dataset()
    })

    observeEvent(input$confirm_convert, {
      removeModal()
      convert_active_dataset()
    })

    convert_active_dataset <- function() {
      ds <- datasets()
      aid <- active_id()
      if (is.null(aid) || is.null(ds[[aid]])) return()
      timed <- tryCatch(
        withProgress(message = "Loading dataset into memory", value = 0.2,
                     with_timing(ensure_in_memory(ds[[aid]]$data))),
        error = function(e) {
          convert_msg(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(timed)) return()
      ds[[aid]]$data <- timed$value
      datasets(ds)
      log_session_event(new_session_event(
        "dataset", "convert_in_memory", dataset = ds[[aid]], parent = ds[[aid]],
        params = list(seconds = round(timed$seconds, 3))), session = session)
      perf_log(perf_log_append(perf_log(), "performance",
                               "convert to in-memory", timed$seconds))
      convert_msg(list(type = "success",
                       message = sprintf("Converted '%s' to in-memory in %.2fs.",
                                         ds[[aid]]$name, timed$seconds)))
    }

    output$convert_status <- renderUI({
      msg <- convert_msg()
      if (is.null(msg)) return(NULL)
      status_alert(msg, extra_class = "mt-2")
    })

    output$timing_table <- renderTable({
      df <- perf_log()
      if (is.null(df) || nrow(df) == 0) {
        return(data.frame(Message = "No operations recorded yet."))
      }
      df <- df[rev(seq_len(nrow(df))), , drop = FALSE]  # most recent first
      names(df) <- c("Time", "Module", "Operation", "Seconds", "Note")
      df
    })
  })
}
