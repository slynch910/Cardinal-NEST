# Cardinal NEST - export module: rds / imzML / Analyze 7.5 / CSV downloads
# Metadata a format cannot embed travels in Cardinal-native .pdata/.fdata/
# .metadata sidecars (writeImzML writes its own; Analyze/CSV are written here).

# Format registry; bundle formats write several files zipped on download.
export_formats <- list(
  rds   = list(label = "R data (.rds)",       ext = "rds",   bundle = FALSE),
  imzML = list(label = "imzML (.imzML/.ibd)", ext = "imzML", bundle = TRUE),
  img   = list(label = "Analyze 7.5 (.img)",  ext = "img",   bundle = TRUE),
  csv   = list(label = "CSV matrix (.csv)",   ext = "csv",   bundle = FALSE)
)

# TRUE when the format writes multiple files that must be zipped for download.
export_is_bundle <- function(format) {
  isTRUE(export_formats[[format]]$bundle)
}

# Legacy sidecar; opt in with options(cardinalnest.legacy_metadata_sidecar).
export_metadata_suffix <- ".metadata.csv"

# x/y/z coordinate columns of pixelData (not user/analysis metadata).
export_coord_names <- function(object) {
  tryCatch(names(as.data.frame(Cardinal::coord(object))),
           error = function(e) character(0))
}

# Added pixelData columns (segment/cluster/ROI/score from Analyze/Spatial).
added_metadata_columns <- function(object) {
  pd <- tryCatch(as.data.frame(Cardinal::pData(object)), error = function(e) NULL)
  if (is.null(pd) || ncol(pd) == 0) return(character(0))
  setdiff(names(pd), c(export_coord_names(object), "run"))
}

# Coord + run + added metadata, one row per pixel; NULL when nothing to
# preserve.
pixel_metadata_frame <- function(object) {
  added <- added_metadata_columns(object)
  if (length(added) == 0) return(NULL)
  pd <- as.data.frame(Cardinal::pData(object))
  ids <- intersect(c(export_coord_names(object), "run"), names(pd))
  pd[, unique(c(ids, added)), drop = FALSE]
}

# Legacy metadata sidecar CSV; NULL when nothing to preserve.
write_metadata_sidecar <- function(object, file) {
  meta <- pixel_metadata_frame(object)
  if (is.null(meta)) return(NULL)
  utils::write.csv(meta, file, row.names = FALSE)
  file
}

# Coords (readAnalyze rebuilds a regular grid, so rows are reconciled by
# coordinate) + run + added columns. Same condition as Cardinal's
# .write_pixelData: nrun > 1 or extra variables exist.
pdata_sidecar_frame <- function(object) {
  pd <- as.data.frame(Cardinal::pData(object))
  coords <- export_coord_names(object)
  vars <- setdiff(names(pd), coords)
  if (!(length(Cardinal::runNames(object)) > 1L || length(vars) > 1L)) {
    return(NULL)
  }
  pd[, unique(c(coords, vars)), drop = FALSE]
}

# mz (restored positionally, undoing binary-format truncation) + added
# featureData columns; NULL when only mz exists.
fdata_sidecar_frame <- function(object) {
  fd <- tryCatch(as.data.frame(Cardinal::fData(object)),
                 error = function(e) NULL)
  if (is.null(fd)) return(NULL)
  vars <- setdiff(names(fd), "mz")
  if (length(vars) == 0) return(NULL)
  fd[, unique(c("mz", vars)), drop = FALSE]
}

# Same write.table/dput output as Cardinal so the files read back identically.
export_native_sidecars <- function(object, file) {
  stem <- tools::file_path_sans_ext(file)
  out <- character(0)

  pd <- pdata_sidecar_frame(object)
  if (!is.null(pd)) {
    p <- paste0(stem, ".pdata")
    utils::write.table(pd, file = p)
    out <- c(out, p)
  }

  fd <- fdata_sidecar_frame(object)
  if (!is.null(fd)) {
    f <- paste0(stem, ".fdata")
    utils::write.table(fd, file = f)
    out <- c(out, f)
  }

  md <- tryCatch(S4Vectors::metadata(object), error = function(e) NULL)
  if (length(md) > 0) {
    m <- paste0(stem, ".metadata")
    ok <- try(dput(md, file = m), silent = TRUE)
    if (!inherits(ok, "try-error")) out <- c(out, m)
  }
  out
}

# Sidecar names an export will produce; none for rds (embedded in the object).
export_sidecar_names <- function(object, base, format) {
  if (is.null(object) || identical(format, "rds")) return(character(0))
  nms <- character(0)
  if (!is.null(pdata_sidecar_frame(object))) {
    nms <- c(nms, paste0(base, ".pdata"))
  }
  if (!is.null(fdata_sidecar_frame(object))) {
    nms <- c(nms, paste0(base, ".fdata"))
  }
  md <- tryCatch(S4Vectors::metadata(object), error = function(e) NULL)
  if (length(md) > 0) nms <- c(nms, paste0(base, ".metadata"))
  if (isTRUE(getOption("cardinalnest.legacy_metadata_sidecar", FALSE)) &&
      length(added_metadata_columns(object)) > 0) {
    nms <- c(nms, paste0(base, export_metadata_suffix))
  }
  nms
}

write_provenance_sidecar <- function(payload, file) {
  if (is.null(payload)) return(NULL)
  writeLines(provenance_json(payload), file)
  file
}

# rds provenance lives in the metadata slot; ISO timestamps keep it identical
# to the JSON sidecar content.
embed_provenance <- function(object, payload) {
  if (is.null(payload)) return(object)
  S4Vectors::metadata(object)$cardinal_nest <- json_safe(payload)
  object
}

# TRUE when the download is a zip: bundle formats, or CSV with sidecars or
# provenance. rds is self-contained and never bundles.
export_delivers_bundle <- function(format, object = NULL, provenance = NULL) {
  if (export_is_bundle(format)) return(TRUE)
  if (identical(format, "csv") &&
      (!is.null(provenance) ||
       (!is.null(object) &&
        length(export_sidecar_names(object, "x", format)) > 0))) {
    return(TRUE)
  }
  FALSE
}

# One row per lineage step, oldest first.
provenance_chain_table <- function(chain) {
  if (length(chain) == 0) {
    return(data.frame(Message = "No provenance recorded for this dataset.",
                      stringsAsFactors = FALSE))
  }
  fmt_time <- function(t) {
    out <- tryCatch(format(t, "%Y-%m-%d %H:%M:%S"), error = function(e) NULL)
    if (is.null(out) || is.na(out)) "" else out
  }
  data.frame(
    Step = seq_along(chain),
    Operation = vapply(chain, function(p) p$operation %||% "", character(1)),
    Method = vapply(chain, function(p) p$method %||% "", character(1)),
    Parent = vapply(chain, function(p) p$parent_name %||% "", character(1)),
    Seed = vapply(chain, function(p) {
      if (is.null(p$seed) || is.na(p$seed)) "" else as.character(p$seed)
    }, character(1)),
    Time = vapply(chain, function(p) fmt_time(p$timestamp), character(1)),
    stringsAsFactors = FALSE
  )
}

# Safe filesystem base name from a user-supplied filename.
export_basename <- function(filename, default = "dataset") {
  nm <- if (is.null(filename)) "" else trimws(filename)
  nm <- tools::file_path_sans_ext(basename(nm))
  nm <- gsub("[^A-Za-z0-9_.-]", "_", nm)
  if (!nzchar(nm)) nm <- default
  nm
}

# Feature-by-pixel CSV with an m/z first column; pixel columns named by
# coordinate when counts match.
export_csv_matrix <- function(object, file) {
  mat <- as.matrix(Cardinal::spectra(object))
  df <- data.frame(mz = Cardinal::mz(object), mat, check.names = FALSE)

  cd <- tryCatch(as.data.frame(Cardinal::coord(object)), error = function(e) NULL)
  if (!is.null(cd) && nrow(cd) == ncol(mat)) {
    px_names <- apply(cd, 1, function(r) paste0("px_", paste(r, collapse = "_")))
    colnames(df)[-1] <- make.unique(px_names)
  }

  utils::write.csv(df, file, row.names = FALSE)
  file
}

# Materialize file-backed spectra so the .rds is self-contained (saveRDS only
# records references to external .ibd/.img files). Column-block reads stay
# under matter's 32-bit size-overflow point.
export_materialize_spectra <- function(object, max_elems = 1e8) {
  s <- tryCatch(Cardinal::spectra(object), error = function(e) NULL)
  if (is.null(s) || is.matrix(s)) return(object)
  nr <- nrow(s)
  nc <- ncol(s)
  if (is.null(nr) || is.null(nc)) {
    Cardinal::spectra(object) <- as.matrix(s)
    return(object)
  }
  out <- matrix(0, nrow = nr, ncol = nc)
  block_ncol <- max(1L, as.integer(floor(max_elems / max(nr, 1))))
  start <- 1L
  while (start <= nc) {
    end <- min(start + block_ncol - 1L, nc)
    out[, start:end] <- as.matrix(s[, start:end, drop = FALSE])
    start <- end + 1L
  }
  Cardinal::spectra(object) <- out
  object
}

# rds (non-matrix spectra) and csv build a full in-memory copy.
export_memory_preflight <- function(object, format) {
  if (is.null(object)) return(NULL)
  needs_copy <- switch(format,
    rds = !is.matrix(tryCatch(Cardinal::spectra(object), error = function(e) NULL)),
    csv = TRUE,
    FALSE)
  if (!isTRUE(needs_copy)) return(NULL)
  memory_preflight(estimate_spectra_bytes(object),
                   action = sprintf("Exporting as %s", export_formats[[format]]$label),
                   alternative = "Consider exporting as imzML or Analyze 7.5 instead.")
}

# Write `object` to `dir` in `format`; created files detected by diffing the
# dir contents before/after the write. Provenance is embedded for rds, else
# written as <base>.provenance.json.
export_dataset <- function(object, dir, base,
                           format = c("rds", "imzML", "img", "csv"),
                           provenance = NULL,
                           legacy_sidecar = getOption(
                             "cardinalnest.legacy_metadata_sidecar", FALSE)) {
  format <- match.arg(format)
  pf <- export_memory_preflight(object, format)
  if (!is.null(pf)) stop_if_memory_blocked(pf)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  before <- list.files(dir, full.names = TRUE)
  prov_file <- file.path(dir, paste0(base, ".provenance.json"))

  write_legacy <- function() {
    if (isTRUE(legacy_sidecar)) {
      write_metadata_sidecar(
        object, file.path(dir, paste0(base, export_metadata_suffix)))
    }
  }

  switch(format,
    rds = {
      # Materialize so the .rds does not depend on the original .ibd/.img.
      obj <- embed_provenance(export_materialize_spectra(object), provenance)
      saveRDS(obj, file.path(dir, paste0(base, ".rds")))
    },
    imzML = {
      # writeImzML writes the .pdata/.fdata/.metadata sidecars itself.
      Cardinal::writeImzML(object, file.path(dir, paste0(base, ".imzML")),
                           bundle = FALSE, verbose = FALSE)
      write_legacy()
      write_provenance_sidecar(provenance, prov_file)
    },
    img = {
      Cardinal::writeAnalyze(object, file.path(dir, paste0(base, ".img")),
                             verbose = FALSE)
      # writeAnalyze writes no sidecars itself.
      export_native_sidecars(object, file.path(dir, paste0(base, ".img")))
      write_legacy()
      write_provenance_sidecar(provenance, prov_file)
    },
    csv = {
      export_csv_matrix(object, file.path(dir, paste0(base, ".csv")))
      export_native_sidecars(object, file.path(dir, paste0(base, ".csv")))
      write_legacy()
      write_provenance_sidecar(provenance, prov_file)
    }
  )

  after <- list.files(dir, full.names = TRUE)
  setdiff(after, before)
}

# Write `format` into a temp dir, then copy a single file or zip multiple
# files into `outfile`.
write_export_download <- function(object, base, format, outfile,
                                  provenance = NULL) {
  work <- tempfile(pattern = "cardinal_export_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  files <- export_dataset(object, work, base, format, provenance = provenance)

  if (export_is_bundle(format) || length(files) > 1) {
    zip::zipr(zipfile = outfile, files = files)
  } else {
    file.copy(files[1], outfile, overwrite = TRUE)
  }
  outfile
}

exportUI <- function(id) {
  ns <- NS(id)

  format_choices <- stats::setNames(names(export_formats),
                                    vapply(export_formats, `[[`, character(1), "label"))

  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
             h3("Export"),
             p("Export the active dataset to disk. imzML and Analyze 7.5 formats produce several files and are delivered as a .zip archive."),
             uiOutput(ns("active_dataset"))
      )
    ),
    fluidRow(
      column(6,
        bslib::card(
          bslib::card_header("Export dataset"),
          bslib::card_body(
            selectInput(ns("format"), "Format", choices = format_choices,
                        selected = "rds"),
            textInput(ns("filename"), "Filename (without extension)",
                      value = "dataset"),
            downloadButton(ns("download"), "Download",
                           class = "btn-success w-100")
          )
        )
      ),
      column(6,
        bslib::card(
          bslib::card_header("Output"),
          bslib::card_body(uiOutput(ns("preview")))
        )
      )
    ),
    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header("Record provenance"),
          bslib::card_body(
            p(class = "text-muted small",
              paste("Processing history of the active dataset, oldest step",
                    "first. This lineage is what the .provenance.json sidecar",
                    "(or the .rds metadata slot) carries on export.")),
            div(style = "overflow-x: auto;",
                tableOutput(ns("provenance_table"))),
            uiOutput(ns("provenance_json"))
          )
        )
      )
    )
  )
}

exportServer <- function(id, datasets, active_id) {
  moduleServer(id, function(input, output, session) {

    mse_data <- active_dataset(datasets, active_id)
    active_record <- active_dataset_record(datasets, active_id)
    output$active_dataset <- active_dataset_ui(datasets, active_id)

    observeEvent(active_record(), {
      rec <- active_record()
      if (is.null(rec)) return()
      updateTextInput(session, "filename",
                      value = export_basename(rec$name, default = "dataset"))
    })

    observe({
      if (is.null(mse_data())) {
        shinyjs::disable("download")
      } else {
        shinyjs::enable("download")
      }
    })

    output$preview <- renderUI({
      fmt <- input$format
      base <- export_basename(input$filename)
      if (is.null(fmt) || is.null(export_formats[[fmt]])) return(NULL)

      data <- mse_data()
      has_meta <- !is.null(data) &&
        (!is.null(pdata_sidecar_frame(data)) || !is.null(fdata_sidecar_frame(data)))

      has_prov <- !is.null(active_record())

      files <- switch(fmt,
        rds   = paste0(base, ".rds"),
        csv   = paste0(base, ".csv"),
        imzML = c(paste0(base, ".imzML"), paste0(base, ".ibd")),
        img   = c(paste0(base, ".hdr"), paste0(base, ".img"), paste0(base, ".t2m"))
      )
      # The binary writers leave a .log next to the data files.
      if (fmt %in% c("imzML", "img")) {
        files <- c(files, paste0(base, ".log"))
      }
      files <- c(files, export_sidecar_names(data, base, fmt))
      if (has_prov && fmt != "rds") {
        files <- c(files, paste0(base, ".provenance.json"))
      }

      delivered <- if (export_delivers_bundle(fmt, data, has_prov)) paste0(base, "_", fmt, ".zip")
                   else files

      meta_note <- if (fmt == "rds" && has_meta) {
        helpText("Pixel and feature metadata is stored inside the .rds object.")
      } else if (has_meta) {
        sidecars <- export_sidecar_names(data, base, fmt)
        helpText(sprintf("Metadata that %s cannot embed is written to %s and restored on re-import.",
                         fmt, paste(sidecars, collapse = ", ")))
      }

      prov_note <- if (has_prov && fmt == "rds") {
        helpText("Record provenance is embedded in the .rds object's metadata slot.")
      } else if (has_prov) {
        helpText(sprintf("Record provenance is written to %s.",
                         paste0(base, ".provenance.json")))
      }

      pf <- export_memory_preflight(data, fmt)
      mem_note <- if (!is.null(pf) && pf$level %in% c("warn", "block")) {
        status_alert(list(type = if (pf$level == "block") "danger" else "warning",
                          message = pf$message), extra_class = "mt-2")
      }

      tagList(
        strong("Files written:"),
        tags$ul(lapply(files, tags$li)),
        strong("Downloaded as: "), code(delivered),
        meta_note,
        prov_note,
        mem_note
      )
    })

    output$provenance_table <- renderTable({
      if (is.null(active_record())) {
        return(data.frame(Message = "No active dataset."))
      }
      provenance_chain_table(record_lineage(datasets(), active_id()))
    }, rownames = FALSE)

    output$provenance_json <- renderUI({
      if (is.null(active_record())) return(NULL)
      payload <- export_provenance(datasets(), active_id(), input$format)
      if (is.null(payload)) return(NULL)
      tags$details(
        tags$summary("Full export payload (JSON)"),
        tags$pre(class = "small mt-2",
                 style = "max-height: 320px; overflow: auto;",
                 provenance_json(payload))
      )
    })

    download_filename <- function() {
      fmt <- input$format
      base <- export_basename(input$filename)
      has_prov <- !is.null(active_record())
      if (export_delivers_bundle(fmt, mse_data(), has_prov)) {
        paste0(base, "_", fmt, ".zip")
      } else {
        paste0(base, ".", export_formats[[fmt]]$ext)
      }
    }

    output$download <- downloadHandler(
      filename = download_filename,
      content = function(file) {
        data <- mse_data()
        if (is.null(data)) {
          stop("No active dataset to export.")
        }
        prov <- export_provenance(datasets(), active_id(), input$format)
        tryCatch(
          write_export_download(data, export_basename(input$filename),
                                input$format, file, provenance = prov),
          cardinalnest_memory_blocked = function(e) {
            showNotification(conditionMessage(e), type = "error", duration = 15)
            stop(e)
          }
        )
        rec <- active_record()
        log_export_event("export_dataset", download_filename(), input$format,
                         dataset = rec, provenance = rec$provenance,
                         session = session)
      }
    )
  })
}
