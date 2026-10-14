# Cardinal NEST - export module: rds / imzML / Analyze 7.5 / CSV downloads

library(shiny)

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

# Suffix used for the per-pixel metadata sidecar written alongside formats that
# cannot embed arbitrary pixelData columns themselves.
export_metadata_suffix <- ".metadata.csv"

# x/y/z coordinate columns of pixelData (not user/analysis metadata).
export_coord_names <- function(object) {
  tryCatch(names(as.data.frame(Cardinal::coord(object))),
           error = function(e) character(0))
}

# Added pixelData columns (segment/cluster/ROI/score from Analyze/Spatial);
# binary formats cannot embed them, so they ship in a sidecar CSV.
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

# Write the metadata sidecar CSV; NULL when nothing to preserve.
write_metadata_sidecar <- function(object, file) {
  meta <- pixel_metadata_frame(object)
  if (is.null(meta)) return(NULL)
  utils::write.csv(meta, file, row.names = FALSE)
  file
}

# TRUE when the download is a multi-file zip (bundle formats, or CSV with a
# metadata sidecar). rds is self-contained and never bundles.
export_delivers_bundle <- function(format, object = NULL) {
  if (export_is_bundle(format)) return(TRUE)
  if (identical(format, "csv") && !is.null(object) &&
      length(added_metadata_columns(object)) > 0) {
    return(TRUE)
  }
  FALSE
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

# Write `object` to `dir` in `format`; created files detected by diffing the
# dir contents before/after the write.
export_dataset <- function(object, dir, base,
                           format = c("rds", "imzML", "img", "csv")) {
  format <- match.arg(format)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  before <- list.files(dir, full.names = TRUE)

  switch(format,
    rds = {
      # rds serializes the whole object including pixelData; materialize
      # first so it does not depend on the original .ibd/.img files.
      saveRDS(export_materialize_spectra(object),
              file.path(dir, paste0(base, ".rds")))
    },
    imzML = {
      Cardinal::writeImzML(object, file.path(dir, paste0(base, ".imzML")),
                           bundle = FALSE, verbose = FALSE)
      write_metadata_sidecar(object,
                             file.path(dir, paste0(base, export_metadata_suffix)))
    },
    img = {
      Cardinal::writeAnalyze(object, file.path(dir, paste0(base, ".img")),
                             verbose = FALSE)
      write_metadata_sidecar(object,
                             file.path(dir, paste0(base, export_metadata_suffix)))
    },
    csv = {
      export_csv_matrix(object, file.path(dir, paste0(base, ".csv")))
      write_metadata_sidecar(object,
                             file.path(dir, paste0(base, export_metadata_suffix)))
    }
  )

  after <- list.files(dir, full.names = TRUE)
  setdiff(after, before)
}

# Write `format` into a temp dir, then copy a single file or zip multiple
# files into `outfile`.
write_export_download <- function(object, base, format, outfile) {
  work <- tempfile(pattern = "cardinal_export_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  files <- export_dataset(object, work, base, format)

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
    )
  )
}

exportServer <- function(id, datasets, active_id) {
  moduleServer(id, function(input, output, session) {
    library(Cardinal)

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
      meta_cols <- if (is.null(data)) character(0) else added_metadata_columns(data)
      has_meta <- length(meta_cols) > 0

      files <- switch(fmt,
        rds   = paste0(base, ".rds"),
        csv   = paste0(base, ".csv"),
        imzML = c(paste0(base, ".imzML"), paste0(base, ".ibd")),
        img   = c(paste0(base, ".hdr"), paste0(base, ".img"), paste0(base, ".t2m"))
      )
      # Binary formats can't embed added pixelData; it ships as a sidecar CSV.
      if (has_meta && fmt != "rds") {
        files <- c(files, paste0(base, export_metadata_suffix))
      }

      delivered <- if (export_delivers_bundle(fmt, data)) paste0(base, "_", fmt, ".zip")
                   else files

      meta_note <- if (fmt == "rds" && has_meta) {
        helpText(sprintf("Added metadata columns (%s) are stored inside the .rds object.",
                         paste(meta_cols, collapse = ", ")))
      } else if (has_meta) {
        helpText(sprintf("Added metadata columns (%s) are written to %s.",
                         paste(meta_cols, collapse = ", "),
                         paste0(base, export_metadata_suffix)))
      }

      tagList(
        strong("Files written:"),
        tags$ul(lapply(files, tags$li)),
        strong("Downloaded as: "), code(delivered),
        meta_note
      )
    })

    output$download <- downloadHandler(
      filename = function() {
        fmt <- input$format
        base <- export_basename(input$filename)
        if (export_delivers_bundle(fmt, mse_data())) {
          paste0(base, "_", fmt, ".zip")
        } else {
          paste0(base, ".", export_formats[[fmt]]$ext)
        }
      },
      content = function(file) {
        data <- mse_data()
        if (is.null(data)) {
          stop("No active dataset to export.")
        }
        write_export_download(data, export_basename(input$filename),
                              input$format, file)
      }
    )
  })
}
