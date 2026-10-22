# Cardinal NEST - shared helper functions

`%||%` <- function(a, b) if (is.null(a)) b else a

# Comma-separated numbers; with integer = TRUE, tokens may also be inclusive
# a:b ranges. Errors carry a user-facing message prefixed by `label`.
parse_numeric_list <- function(str, label = "Input", integer = FALSE) {
  hint <- if (integer) {
    "Enter whole numbers or a:b ranges separated by commas, e.g. 1:10, 30."
  } else {
    "Enter one or more numbers separated by commas."
  }
  fail <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
  tokens <- trimws(strsplit(paste(str %||% "", collapse = ","), ",", fixed = TRUE)[[1]])
  if (length(tokens) == 0 || all(tokens == "")) fail("%s is empty. %s", label, hint)

  vals <- lapply(tokens, function(tok) {
    if (grepl(":", tok, fixed = TRUE)) {
      if (!integer) fail("%s does not support ranges ('%s'). %s", label, tok, hint)
      ends <- suppressWarnings(as.numeric(trimws(strsplit(tok, ":", fixed = TRUE)[[1]])))
      if (length(ends) != 2 || any(!is.finite(ends)))
        fail("%s has an invalid range '%s'. %s", label, tok, hint)
      if (any(ends < 1 | ends != round(ends)))
        fail("%s must contain positive whole numbers ('%s'). %s", label, tok, hint)
      if (ends[1] > ends[2])
        fail("%s has a reversed range '%s'; use %s:%s instead.",
             label, tok, format(ends[2]), format(ends[1]))
      return(seq(ends[1], ends[2]))
    }
    v <- suppressWarnings(as.numeric(tok))
    if (!is.finite(v)) fail("%s has an invalid entry '%s'. %s", label, tok, hint)
    if (integer && (v < 1 || v != round(v)))
      fail("%s must contain positive whole numbers ('%s'). %s", label, tok, hint)
    v
  })
  unique(as.numeric(unlist(vals)))
}

# Shiny progress wrapper: sit at `mid` while `expr` (a blocking call) runs.
with_progress_stages <- function(message, expr, start = 0.15, mid = 0.55) {
  shiny::withProgress(message = message, value = start, {
    shiny::setProgress(mid, detail = "Computing")
    out <- force(expr)
    shiny::setProgress(1, detail = "Done")
    out
  })
}

# Human-readable byte size, e.g. 2.1e9 -> "2 GB".
format_bytes <- function(bytes) {
  bytes <- as.numeric(bytes)
  if (length(bytes) == 0 || is.na(bytes)) return("NA")
  units <- c("B", "KB", "MB", "GB", "TB", "PB")
  if (bytes < 1) return("0 B")
  power <- min(floor(log(bytes, 1024)), length(units) - 1)
  value <- bytes / (1024^power)
  sprintf("%s %s", formatC(value, format = "fg", digits = 3), units[power + 1])
}

# TRUE when spectra are stored out-of-memory (file-backed/matter).
is_file_backed <- function(data) {
  sp <- tryCatch(Cardinal::spectra(data), error = function(e) NULL)
  !is.null(sp) && !is.matrix(sp)
}

# Cheap O(1) cache key from dimensions and m/z extent only.
dataset_fingerprint <- function(data) {
  if (is.null(data)) return("null")
  mzr <- tryCatch(range(Cardinal::mz(data), na.rm = TRUE),
                  error = function(e) c(NA_real_, NA_real_))
  paste(nrow(data), ncol(data),
        signif(mzr[1], 8), signif(mzr[2], 8), sep = ":")
}

# Wall-clock timing wrapper returning list(value, seconds).
with_timing <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

# In-memory key/value cache backed by an environment.
new_result_cache <- function() {
  store <- new.env(parent = emptyenv())
  list(
    get = function(key) if (exists(key, envir = store, inherits = FALSE)) {
      get(key, envir = store, inherits = FALSE)
    } else NULL,
    set = function(key, value) assign(key, value, envir = store),
    has = function(key) exists(key, envir = store, inherits = FALSE),
    clear = function() rm(list = ls(store, all.names = TRUE), envir = store),
    size = function() length(ls(store, all.names = TRUE))
  )
}

# Active dataset's MSImagingExperiment object (or NULL).
active_dataset <- function(datasets, active_id) {
  shiny::reactive({
    ds <- datasets()
    aid <- active_id()
    if (is.null(aid) || is.null(ds[[aid]])) return(NULL)
    ds[[aid]]$data
  })
}

# Full active dataset record (or NULL).
active_dataset_record <- function(datasets, active_id) {
  shiny::reactive({
    ds <- datasets()
    aid <- active_id()
    if (is.null(aid) || is.null(ds[[aid]])) return(NULL)
    ds[[aid]]
  })
}

# "Active dataset: <name>" indicator; warning alert when none is active.
active_dataset_ui <- function(datasets, active_id) {
  shiny::renderUI({
    ds <- datasets()
    aid <- active_id()
    if (is.null(aid) || is.null(ds[[aid]])) {
      return(shiny::div(class = "alert alert-warning",
        "No active dataset. Import or select a dataset in the Import tab."))
    }
    shiny::div(class = "alert alert-info",
               shiny::strong("Active dataset: "), ds[[aid]]$name)
  })
}

# Bootstrap alert class for a status type.
status_alert_class <- function(type) {
  switch(type %||% "info",
         success = "alert-success",
         danger = "alert-danger",
         warning = "alert-warning",
         "alert-info")
}

# Bootstrap alert div from a status list(type, message).
status_alert <- function(st, extra_class = NULL) {
  cls <- trimws(paste("alert", status_alert_class(st$type), extra_class))
  shiny::div(class = cls, st$message)
}

# renderUI wrapper for a status reactive.
render_status <- function(status_reactive) {
  shiny::renderUI(status_alert(status_reactive()))
}

# Standard "save as new dataset" modal; `confirm_id` is the (un-namespaced)
# confirm-button id, filename_id/ext_id overridable for modules with several.
save_dataset_modal <- function(ns, title, default_name, confirm_id,
                               confirm_label = "Save",
                               filename_id = "save_filename",
                               ext_id = "save_ext") {
  shiny::showModal(shiny::modalDialog(
    title = title,
    shiny::textInput(ns(filename_id), "Filename", value = default_name),
    shiny::selectInput(ns(ext_id), "Extension",
      choices = c("R data (.rds)" = "rds",
                  "imzML (.imzML)" = "imzML",
                  "Analyze 7.5 (.img)" = "img",
                  "CSV (.csv)" = "csv")),
    footer = shiny::tagList(
      shiny::modalButton("Cancel"),
      shiny::actionButton(ns(confirm_id), confirm_label, class = "btn-primary")
    )
  ))
}

# downloadHandler exporting res$plot_fn() as PDF or PNG; `basename` may be a
# string or a zero-arg function.
export_plot_handler <- function(result_reactive, format_fn, basename,
                                empty_msg = "No result to export.",
                                width = 8, height = 6) {
  shiny::downloadHandler(
    filename = function() {
      stem <- if (is.function(basename)) basename() else basename
      paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", format_fn())
    },
    content = function(file) {
      res <- result_reactive()
      if (identical(format_fn(), "pdf")) {
        grDevices::pdf(file, width = width, height = height)
      } else {
        grDevices::png(file, width = 1600, height = 1200, res = 200)
      }
      on.exit(grDevices::dev.off(), add = TRUE)
      if (is.null(res)) {
        graphics::plot.new(); graphics::text(0.5, 0.5, empty_msg)
      } else {
        res$plot_fn()
      }
    }
  )
}

# Attach per-pixel columns to `data` and build a dataset record. Shared by
# the analysis and spatial save flows.
make_result_record <- function(data, columns, id, label,
                               filename = NULL, ext = NULL) {
  d <- data
  for (nm in names(columns)) {
    Cardinal::pData(d)[[nm]] <- columns[[nm]]
  }
  make_preproc_record(d, id, label, filename = filename, ext = ext)
}

# Append a timing row to a perf-log data.frame.
perf_log_append <- function(df, module, operation, seconds, note = "") {
  row <- data.frame(
    time = format(Sys.time(), "%H:%M:%S"),
    module = module,
    operation = operation,
    seconds = round(as.numeric(seconds), 3),
    note = note,
    stringsAsFactors = FALSE
  )
  if (is.null(df) || nrow(df) == 0) row else rbind(df, row)
}
