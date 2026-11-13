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

# Estimated in-memory spectra size: n_features * n_pixels * 8 bytes.
estimate_spectra_bytes <- function(data) {
  as.numeric(nrow(data)) * as.numeric(ncol(data)) * 8
}

# Bytes, or a K/M/G/T (1024-based) suffix like "512M"; NA when invalid.
parse_memory_size <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return(NA_real_)
  num <- suppressWarnings(as.numeric(x))
  if (is.na(num)) {
    s <- toupper(gsub("\\s+", "", as.character(x)))
    if (!grepl("^[0-9]*\\.?[0-9]+[KMGT]?B?$", s)) return(NA_real_)
    unit <- sub("^[0-9.]+([KMGT]?)B?$", "\\1", s)
    num <- as.numeric(sub("^([0-9.]+).*$", "\\1", s)) *
      1024^(match(unit, c("", "K", "M", "G", "T")) - 1L)
  }
  if (!is.finite(num) || num < 0) NA_real_ else num
}

system_available_memory <- function() {
  tryCatch(as.numeric(ps::ps_system_memory()[["avail"]]),
           error = function(e) NA_real_)
}

session_rss_bytes <- function() {
  tryCatch(as.numeric(ps::ps_memory_info()[["rss"]]),
           error = function(e) NA_real_)
}

# Tightest cgroup (v1 or v2) limit-minus-usage across every level of this
# process's cgroup path; NA when no level sets a limit. Containers and
# cluster schedulers enforce memory this way.
cgroup_memory_headroom <- function(root = "/sys/fs/cgroup",
                                   proc_cgroup = "/proc/self/cgroup") {
  read_num <- function(f) {
    if (!file.exists(f)) return(NA_real_)
    v <- suppressWarnings(as.numeric(readLines(f, n = 1L, warn = FALSE)))
    if (length(v) == 1L && is.finite(v)) v else NA_real_
  }
  levels_of <- function(path) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    c(vapply(rev(seq_along(parts)),
             function(i) paste(parts[seq_len(i)], collapse = "/"), character(1)), "")
  }
  headroom <- function(dirs, limit_file, usage_file) {
    vals <- vapply(dirs, function(d) {
      lim <- read_num(file.path(d, limit_file))
      if (is.na(lim) || lim >= 2^60) return(NA_real_)
      use <- read_num(file.path(d, usage_file))
      max(lim - (if (is.na(use)) 0 else use), 0)
    }, numeric(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) NA_real_ else min(vals)
  }
  lines <- if (file.exists(proc_cgroup)) readLines(proc_cgroup, warn = FALSE) else character(0)
  v2 <- sub("^0::", "", grep("^0::", lines, value = TRUE))
  v1_pattern <- "^[0-9]+:([^:]*,)?memory(,[^:]*)?:"
  v1 <- sub(v1_pattern, "", grep(v1_pattern, lines, value = TRUE))
  out <- c(
    headroom(file.path(root, levels_of(if (length(v2)) v2[1] else "")),
             "memory.max", "memory.current"),
    headroom(file.path(root, "memory", levels_of(if (length(v1)) v1[1] else "")),
             "memory.limit_in_bytes", "memory.usage_in_bytes")
  )
  out <- out[!is.na(out)]
  if (length(out) == 0) NA_real_ else min(out)
}

# Slurm job memory (bytes): SLURM_MEM_PER_NODE, else per-CPU x CPUs (MB).
slurm_memory_limit <- function(env = Sys.getenv()) {
  mb <- function(nm) suppressWarnings(as.numeric(unname(env[nm])))
  node <- mb("SLURM_MEM_PER_NODE")
  if (is.finite(node)) return(node * 1024^2)
  cpu <- mb("SLURM_MEM_PER_CPU")
  if (!is.finite(cpu)) return(NA_real_)
  n <- mb("SLURM_CPUS_ON_NODE")
  cpu * 1024^2 * (if (is.finite(n)) n else 1)
}

# options(cardinalnest.memory_limit) or CARDINALNEST_MEMORY_LIMIT.
user_memory_limit <- function() {
  parse_memory_size(getOption("cardinalnest.memory_limit",
                              Sys.getenv("CARDINALNEST_MEMORY_LIMIT", "")))
}

# Smallest of system available, cgroup headroom, and Slurm/user limits
# minus current session RSS.
available_memory_bytes <- function(system = system_available_memory(),
                                   cgroup = cgroup_memory_headroom(),
                                   slurm = slurm_memory_limit(),
                                   limit = user_memory_limit(),
                                   rss = session_rss_bytes()) {
  used <- if (is.na(rss)) 0 else rss
  totals <- c(slurm, limit)
  totals <- totals[!is.na(totals)]
  vals <- c(system, cgroup, pmax(totals - used, 0))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) NA_real_ else min(vals)
}

# level: "ok" | "warn" (> warn_frac of available) | "block" (> block_frac)
# | "unknown". `alternative` is appended to warn/block messages.
memory_preflight <- function(bytes, available = available_memory_bytes(),
                             warn_frac = getOption("cardinalnest.memory_warn_frac", 0.5),
                             block_frac = getOption("cardinalnest.memory_block_frac", 0.8),
                             action = "Loading the spectra into memory",
                             alternative = NULL) {
  bytes <- as.numeric(bytes)
  available <- as.numeric(available)
  if (length(available) != 1L || is.na(available)) {
    return(list(level = "unknown", bytes = bytes, available = NA_real_,
                message = sprintf("%s needs about %s; available memory could not be determined.",
                                  action, format_bytes(bytes))))
  }
  frac <- if (bytes <= 0) 0 else if (available <= 0) Inf else bytes / available
  level <- if (frac > block_frac) "block" else if (frac > warn_frac) "warn" else "ok"
  msg <- sprintf("%s needs about %s of the %s available%s.", action,
                 format_bytes(bytes), format_bytes(available),
                 if (is.finite(frac)) sprintf(" (%.0f%%)", 100 * frac) else "")
  msg <- switch(level,
    block = paste(msg, sprintf("Blocked: this exceeds %.0f%% of available memory.",
                               100 * block_frac), alternative),
    warn = paste(msg, sprintf("This exceeds %.0f%% of available memory and may slow down or crash the session.",
                              100 * warn_frac), alternative),
    msg)
  list(level = level, bytes = bytes, available = available, message = trimws(msg))
}

stop_if_memory_blocked <- function(preflight) {
  if (identical(preflight$level, "block")) {
    stop(structure(class = c("cardinalnest_memory_blocked", "error", "condition"),
                   list(message = preflight$message, call = NULL,
                        preflight = preflight)))
  }
  invisible(preflight)
}

# Stop on "block"; call on_warn(message) on "warn".
memory_guard <- function(bytes, ..., on_warn = NULL) {
  pf <- stop_if_memory_blocked(memory_preflight(bytes, ...))
  if (identical(pf$level, "warn") && is.function(on_warn)) on_warn(pf$message)
  invisible(pf)
}

notify_memory_warning <- function(message) {
  shiny::showNotification(message, type = "warning", duration = 15)
}

# "file" | "shared" (matter path "@<uuid>") | "memory"; sparse matrices are
# resolved through their data slot.
spectra_storage <- function(data) {
  sp <- tryCatch(Cardinal::spectra(data), error = function(e) NULL)
  if (is.null(sp)) return(NA_character_)
  if (methods::is(sp, "sparse_mat")) sp <- sp@data
  if (!methods::is(sp, "matter")) return("memory")
  p <- tryCatch(matter::path(sp), error = function(e) NULL)
  if (length(p) == 0 || all(is.na(p))) return("memory")
  if (all(startsWith(p, "@"))) "shared" else "file"
}

is_file_backed <- function(data) {
  identical(spectra_storage(data), "file")
}

storage_label <- function(data) {
  switch(spectra_storage(data) %||% "unknown",
         file = "File-backed",
         shared = "Shared memory",
         memory = "In-memory",
         "Unknown")
}

# Store for file-backed derived datasets; session tempdir unless
# options(cardinalnest.store_dir) is set.
nest_store_dir <- function() {
  dir <- getOption("cardinalnest.store_dir",
                   file.path(tempdir(), "cardinal_nest_store"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  normalizePath(dir, winslash = "/")
}

nest_store_file <- function(prefix = "dataset") {
  prefix <- gsub("[^A-Za-z0-9_-]", "_", prefix)
  tempfile(pattern = paste0(prefix, "_"), tmpdir = nest_store_dir(),
           fileext = ".imzML")
}

# Chunked summarizeFeatures() sum propagates NA without loading the matrix.
spectra_has_na <- function(data) {
  sp <- tryCatch(Cardinal::spectra(data), error = function(e) NULL)
  if (is.null(sp)) return(FALSE)
  if (is.matrix(sp)) return(anyNA(sp))
  s <- Cardinal::summarizeFeatures(data, stat = "sum", verbose = FALSE)
  anyNA(Cardinal::fData(s)[["sum"]])
}

# Write to the store as imzML and re-read file-backed, keeping metadata.
persist_file_backed <- function(data, prefix = "dataset") {
  file <- nest_store_file(prefix)
  Cardinal::writeMSIData(data, file)
  out <- Cardinal::readMSIData(file, memory = FALSE, verbose = FALSE)
  if (!methods::is(out, "MSImagingExperiment")) {
    out <- methods::as(out, "MSImagingExperiment")
  }
  Cardinal::pixelData(out) <- Cardinal::pixelData(data)
  Cardinal::featureData(out) <- Cardinal::featureData(data)
  out
}

# Persist `result` file-backed when `source` was.
keep_file_backed <- function(result, source, prefix = "dataset") {
  if (is_file_backed(source) && !is_file_backed(result)) {
    persist_file_backed(result, prefix)
  } else {
    result
  }
}

# matter's chunked apply needs L'Ecuyer-CMRG for parallel-safe streams.
nest_rng_kind <- "L'Ecuyer-CMRG"

# Integer in [0, .Machine$integer.max], else NA.
parse_seed <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return(NA_integer_)
  v <- suppressWarnings(as.numeric(trimws(as.character(x))))
  if (!is.finite(v) || v < 0 || v > .Machine$integer.max || v != round(v)) {
    return(NA_integer_)
  }
  as.integer(v)
}

# options(cardinalnest.seed) or CARDINALNEST_SEED, else 42.
default_seed <- function() {
  s <- parse_seed(getOption("cardinalnest.seed", Sys.getenv("CARDINALNEST_SEED", "")))
  if (is.na(s)) 42L else s
}

generate_seed <- function() {
  sample.int(.Machine$integer.max, 1L)
}

# `value`, else `global`, else default_seed().
resolve_seed <- function(value, global = NULL) {
  s <- parse_seed(value)
  if (is.na(s)) s <- parse_seed(global)
  if (is.na(s)) default_seed() else s
}

session_global_seed <- function(perf_state = NULL) {
  g <- if (is.null(perf_state)) NULL else perf_state$global_seed
  resolve_seed(g)
}

# Seeded evaluation restoring the caller's RNG state; NULL seed = unseeded.
with_nest_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  withr::with_seed(as.integer(seed), expr, .rng_kind = nest_rng_kind)
}

seed_label <- function(seed) {
  if (is.null(seed)) "" else sprintf(" (seed %d)", as.integer(seed))
}

seed_tag <- function(seed) {
  if (is.null(seed)) "" else sprintf("_seed%d", as.integer(seed))
}

# Seed input plus a "generate" button; wire with observe_seed_generators().
seed_input <- function(ns, id, label = "Random seed (blank = global seed)") {
  shiny::div(class = "d-flex gap-2 align-items-end",
    shiny::div(class = "flex-grow-1",
      shiny::numericInput(ns(id), label, value = NA, min = 0, step = 1)),
    shiny::actionButton(ns(paste0(id, "_generate")), NULL,
                        icon = shiny::icon("dice"),
                        title = "Generate a new seed",
                        class = "btn-outline-secondary mb-3")
  )
}

observe_seed_generators <- function(input, session, ids) {
  lapply(ids, function(id) {
    shiny::observeEvent(input[[paste0(id, "_generate")]], {
      shiny::updateNumericInput(session, id, value = generate_seed())
    }, ignoreInit = TRUE)
  })
  invisible(NULL)
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
                                width = 8, height = 6,
                                session = shiny::getDefaultReactiveDomain()) {
  force(session)
  plot_filename <- function() {
    stem <- if (is.function(basename)) basename() else basename
    paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", format_fn())
  }
  shiny::downloadHandler(
    filename = plot_filename,
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
        log_export_event("export_plot", plot_filename(), format_fn(),
                         provenance = res$provenance, session = session)
      }
    }
  )
}

# Attach per-pixel columns to `data` and build a dataset record. Shared by
# the analysis and spatial save flows.
make_result_record <- function(data, columns, id, label,
                               filename = NULL, ext = NULL,
                               provenance = NULL) {
  d <- data
  for (nm in names(columns)) {
    if (length(columns[[nm]]) != ncol(d)) {
      stop(sprintf("Result column '%s' has %d values but the dataset has %d pixels.",
                   nm, length(columns[[nm]]), ncol(d)), call. = FALSE)
    }
    Cardinal::pData(d)[[nm]] <- columns[[nm]]
  }
  make_preproc_record(d, id, label, filename = filename, ext = ext,
                      provenance = provenance)
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

# Installed package version, else DESCRIPTION when run from source.
nest_app_version <- function(desc_file = "DESCRIPTION") {
  v <- suppressWarnings(
    tryCatch(as.character(utils::packageVersion("CardinalNEST")),
             error = function(e) NA_character_)
  )
  if (!is.na(v)) return(v)
  if (file.exists(desc_file)) {
    dcf <- tryCatch(read.dcf(desc_file), error = function(e) NULL)
    if (!is.null(dcf) && "Version" %in% colnames(dcf)) {
      return(unname(dcf[1, "Version"]))
    }
  }
  NA_character_
}

nest_versions <- function() {
  ver <- function(pkg) {
    suppressWarnings(
      tryCatch(as.character(utils::packageVersion(pkg)),
               error = function(e) NA_character_)
    )
  }
  list(
    R = R.version.string,
    Cardinal = ver("Cardinal"),
    matter = ver("matter"),
    BiocParallel = ver("BiocParallel"),
    CardinalNEST = nest_app_version()
  )
}

# Cardinal's registered BPPARAM, else BiocParallel's default.
nest_backend_desc <- function(param = NULL) {
  if (is.null(param)) {
    param <- tryCatch(Cardinal::getCardinalBPPARAM(), error = function(e) NULL)
  }
  if (is.null(param)) {
    param <- tryCatch(BiocParallel::bpparam(), error = function(e) NULL)
  }
  if (is.null(param)) {
    return(list(type = "SerialParam", workers = 1L))
  }
  workers <- tryCatch(as.integer(BiocParallel::bpnworkers(param)),
                      error = function(e) NA_integer_)
  list(type = class(param)[1], workers = workers)
}

# `parent` is NULL for imports; `params` must be JSON-serializable.
new_provenance <- function(operation, parent = NULL, method = NULL,
                           params = list(), seed = NULL) {
  list(
    parent_id = if (is.null(parent)) NULL else parent$id,
    parent_name = if (is.null(parent)) NULL else parent$name,
    operation = operation,
    method = method,
    params = params,
    seed = if (is.null(seed) || is.na(seed)) NULL else as.integer(seed),
    backend = nest_backend_desc(),
    versions = nest_versions(),
    timestamp = Sys.time()
  )
}

# Provenance chain oldest-first; stops at missing parents or records
# without provenance.
record_lineage <- function(store, id) {
  out <- list()
  seen <- character(0)
  cur_id <- id
  while (!is.null(cur_id) && !(cur_id %in% seen)) {
    seen <- c(seen, cur_id)
    rec <- store[[cur_id]]
    if (is.null(rec)) break
    prov <- rec$provenance
    if (!is.null(prov)) out <- c(list(prov), out)
    cur_id <- if (!is.null(prov)) prov$parent_id else NULL
  }
  out
}

# POSIXt -> ISO-8601, factor -> character, recursively.
json_safe <- function(x) {
  if (inherits(x, "POSIXt")) return(format(x, "%Y-%m-%dT%H:%M:%S%z"))
  if (is.factor(x)) return(as.character(x))
  if (is.list(x)) return(lapply(x, json_safe))
  x
}

export_provenance <- function(store, id, format = NULL) {
  rec <- store[[id]]
  if (is.null(rec)) return(NULL)
  list(
    kind = "cardinalnest_provenance",
    version = 1L,
    dataset = list(id = rec$id, name = rec$name, format = rec$format),
    export = list(operation = "export",
                  format = format,
                  timestamp = Sys.time(),
                  versions = nest_versions()),
    provenance = rec$provenance,
    lineage = record_lineage(store, id)
  )
}

provenance_json <- function(payload, pretty = TRUE) {
  if (is.null(payload)) return(NULL)
  as.character(jsonlite::toJSON(json_safe(payload), auto_unbox = TRUE,
                                pretty = pretty))
}

# Session operation log: a list of events held in session$userData$op_log.
# An event is a provenance record plus category, dataset, seq, and time;
# a supplied `provenance` is used as-is (`operation` overrides it).
new_session_event <- function(category, operation = NULL, dataset = NULL,
                              parent = NULL, method = NULL, params = list(),
                              seed = NULL, provenance = NULL, note = NULL) {
  prov <- provenance %||%
    new_provenance(operation %||% category, parent = parent, method = method,
                   params = params, seed = seed)
  if (!is.null(operation)) prov$operation <- operation
  c(list(seq = NA_integer_,
         logged_at = Sys.time(),
         category = category,
         dataset_id = dataset$id,
         dataset_name = dataset$name),
    prov,
    list(note = note))
}

oplog_append <- function(log, event) {
  log <- log %||% list()
  event$seq <- length(log) + 1L
  c(log, list(event))
}

# No-op when the session publishes no op_log.
log_session_event <- function(event, session = shiny::getDefaultReactiveDomain()) {
  rv <- if (is.null(session)) NULL else session$userData$op_log
  if (!is.function(rv)) return(invisible(NULL))
  shiny::isolate(rv(oplog_append(rv(), event)))
  invisible(event)
}

log_result_event <- function(result, store, category = "analysis",
                             session = shiny::getDefaultReactiveDomain()) {
  prov <- result$provenance
  if (is.null(prov)) return(invisible(NULL))
  dataset <- if (is.null(prov$parent_id)) NULL else store[[prov$parent_id]]
  dataset <- dataset %||% list(id = prov$parent_id, name = prov$parent_name)
  log_session_event(new_session_event(category, dataset = dataset,
                                      provenance = prov),
                    session = session)
}

# Logs each new result once; view changes keep the provenance, so skip.
observe_result_log <- function(result_rv, datasets, category = "analysis",
                               session = shiny::getDefaultReactiveDomain()) {
  last <- NULL
  shiny::observeEvent(result_rv(), {
    prov <- result_rv()$provenance
    if (is.null(prov) || identical(prov, last)) return()
    last <<- prov
    log_result_event(result_rv(), datasets(), category, session = session)
  }, domain = session)
}

log_export_event <- function(operation, file, format = NULL, dataset = NULL,
                             provenance = NULL, params = list(),
                             session = shiny::getDefaultReactiveDomain()) {
  params <- Filter(Negate(is.null),
                   c(list(file = basename(file), format = format,
                          source_operation = provenance$operation), params))
  if (is.null(dataset) && !is.null(provenance$parent_id)) {
    dataset <- list(id = provenance$parent_id, name = provenance$parent_name)
  }
  log_session_event(new_session_event(
    "export", operation, dataset = dataset,
    parent = dataset, method = format, params = params,
    seed = provenance$seed), session = session)
}

# id/name only, so removed records' data is not retained.
store_snapshot <- function(store) {
  lapply(stats::setNames(names(store), names(store)), function(id) {
    list(id = id, name = store[[id]]$name)
  })
}

# Created/removed/renamed events between a snapshot and the current store.
dataset_store_events <- function(old, new) {
  old <- old %||% list()
  new <- new %||% list()
  events <- list()
  for (id in setdiff(names(new), names(old))) {
    rec <- new[[id]]
    prov <- rec$provenance
    category <- if (identical(prov$operation, "import")) "import" else "dataset"
    events <- c(events, list(new_session_event(
      category,
      operation = if (is.null(prov)) "create_dataset" else NULL,
      dataset = list(id = id, name = rec$name),
      provenance = prov)))
  }
  for (id in setdiff(names(old), names(new))) {
    events <- c(events, list(new_session_event(
      "dataset", "remove_dataset", dataset = old[[id]])))
  }
  for (id in intersect(names(old), names(new))) {
    old_name <- old[[id]]$name
    new_name <- new[[id]]$name
    if (!identical(old_name, new_name)) {
      events <- c(events, list(new_session_event(
        "dataset", "rename_dataset", dataset = list(id = id, name = new_name),
        params = list(old_name = old_name, new_name = new_name))))
    }
  }
  events
}

oplog_params_summary <- function(params, width = 160) {
  if (length(params) == 0) return("")
  txt <- tryCatch(provenance_json(params, pretty = FALSE),
                  error = function(e) "")
  if (nchar(txt) > width) paste0(substr(txt, 1, width - 3), "...") else txt
}

oplog_table <- function(log) {
  if (length(log) == 0) {
    return(data.frame(Message = "No operations recorded in this session yet.",
                      stringsAsFactors = FALSE))
  }
  chr <- function(x) if (is.null(x) || length(x) == 0 || is.na(x[1])) "" else as.character(x[1])
  fmt_time <- function(t) {
    out <- tryCatch(format(t, "%Y-%m-%d %H:%M:%S"), error = function(e) "")
    chr(out)
  }
  data.frame(
    `#` = vapply(log, function(e) as.integer(e$seq %||% NA_integer_), integer(1)),
    Time = vapply(log, function(e) fmt_time(e$logged_at), character(1)),
    Category = vapply(log, function(e) chr(e$category), character(1)),
    Operation = vapply(log, function(e) chr(e$operation), character(1)),
    Method = vapply(log, function(e) chr(e$method), character(1)),
    Dataset = vapply(log, function(e) chr(e$dataset_name), character(1)),
    Parent = vapply(log, function(e) chr(e$parent_name), character(1)),
    Seed = vapply(log, function(e) chr(e$seed), character(1)),
    Parameters = vapply(log, function(e) oplog_params_summary(e$params), character(1)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

oplog_payload <- function(log) {
  list(kind = "cardinalnest_operation_record",
       version = 1L,
       written_at = Sys.time(),
       versions = nest_versions(),
       events = log %||% list())
}

oplog_json <- function(log, pretty = TRUE) {
  provenance_json(oplog_payload(log), pretty = pretty)
}

oplog_text <- function(log) {
  v <- nest_versions()
  header <- c("Cardinal NEST session operation record",
              sprintf("Written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
              sprintf("Versions: %s",
                      paste(names(v), unlist(lapply(v, function(x) x %||% NA)),
                            sep = " ", collapse = "; ")),
              "")
  if (length(log) == 0) return(c(header, "No operations recorded."))
  tbl <- oplog_table(log)
  blank <- function(x) if (nzchar(x)) x else "-"
  body <- unlist(lapply(seq_len(nrow(tbl)), function(i) {
    r <- tbl[i, ]
    e <- log[[i]]
    ds_id <- if (is.null(e$dataset_id)) "" else sprintf(" (%s)", e$dataset_id)
    c(sprintf("[%d] %s  %s  %s%s", r[["#"]], r$Time, r$Category, r$Operation,
              if (nzchar(r$Method)) sprintf(" (%s)", r$Method) else ""),
      sprintf("    dataset: %s%s", blank(r$Dataset), ds_id),
      sprintf("    parent:  %s", blank(r$Parent)),
      sprintf("    seed:    %s", blank(r$Seed)),
      sprintf("    params:  %s", blank(oplog_params_summary(e$params, width = Inf))),
      "")
  }))
  c(header, body)
}
