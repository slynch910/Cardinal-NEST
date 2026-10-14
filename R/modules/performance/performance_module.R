# Cardinal NEST - performance module
# BiocParallel backend selection, dataset memory table, result cache
# controls, and a timing log.

library(shiny)

# MulticoreParam forks, unavailable on Windows.
available_backends <- function(os = .Platform$OS.type) {
  backends <- c("Serial (single core)" = "serial",
                "Snow (multiple processes)" = "snow")
  if (!identical(os, "windows")) {
    backends <- c(backends, "Multicore (fork)" = "multicore")
  }
  backends
}

# BiocParallel param object for the chosen backend and worker count.
make_bpparam <- function(backend = "serial", workers = 1L) {
  backend <- match.arg(backend, c("serial", "snow", "multicore"))
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) workers <- 1L
  switch(backend,
         serial = BiocParallel::SerialParam(),
         snow = BiocParallel::SnowParam(workers = workers),
         multicore = BiocParallel::MulticoreParam(workers = workers))
}

# One-line description of a BiocParallel param object.
describe_bpparam <- function(param) {
  cls <- class(param)[1]
  workers <- tryCatch(BiocParallel::bpnworkers(param),
                      error = function(e) NA_integer_)
  if (is.na(workers)) return(cls)
  sprintf("%s, %d worker%s", cls, workers, if (workers == 1L) "" else "s")
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
      Storage = if (isTRUE(is_file_backed(d))) "File-backed" else "In-memory",
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
                 actionButton(ns("apply_backend"), "Apply backend",
                              icon = icon("microchip"), class = "btn-primary"),
                 uiOutput(ns("backend_status"))
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
             )
      ),
      column(7,
             bslib::card(
               bslib::card_header("Dataset memory"),
               bslib::card_body(
                 tableOutput(ns("memory_table")),
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

    backend_msg <- reactiveVal(NULL)
    convert_msg <- reactiveVal(NULL)

    observeEvent(input$cache_enabled, {
      perf_state$cache_enabled <- isTRUE(input$cache_enabled)
    }, ignoreInit = FALSE)

    observeEvent(input$apply_backend, {
      param <- tryCatch(
        make_bpparam(input$backend, input$workers),
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

    output$convert_ui <- renderUI({
      ds <- datasets()
      aid <- active_id()
      if (is.null(aid) || is.null(ds[[aid]])) return(NULL)
      if (!isTRUE(is_file_backed(ds[[aid]]$data))) {
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
      perf_log(perf_log_append(perf_log(), "performance",
                               "convert to in-memory", timed$seconds))
      convert_msg(list(type = "success",
                       message = sprintf("Converted '%s' to in-memory in %.2fs.",
                                         ds[[aid]]$name, timed$seconds)))
    })

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
