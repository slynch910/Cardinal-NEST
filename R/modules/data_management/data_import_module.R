# Cardinal NEST - data import module: imzML / Analyze 7.5 / CSV / rds

library(shiny)

summarize_msi <- function(data) {
  list(
    n_pixels = ncol(data),
    n_features = nrow(data),
    n_runs = length(Cardinal::runNames(data)),
    mz_range = range(Cardinal::mz(data), na.rm = TRUE),
    centroided = centroided(data),
    class = class(data)[1]
  )
}

# Must match export_metadata_suffix in the export module.
import_metadata_suffix <- ".metadata.csv"

# Regex matching a metadata sidecar file name.
import_metadata_pattern <- "\\.metadata\\.csv$"

# Group files into datasets: one per file (csv/rds) or per shared base name
# (imzML primary = .imzML, analyze primary = .hdr). A *.metadata.csv sidecar
# attaches to its base-name group for pixelData restoration.
import_group_files <- function(file_names, format) {
  is_sidecar <- grepl(import_metadata_pattern, file_names, ignore.case = TRUE)
  sidecars <- file_names[is_sidecar]
  sidecar_bases <- sub(import_metadata_pattern, "", sidecars, ignore.case = TRUE)
  data_files <- file_names[!is_sidecar]

  # Attach the sidecar whose base matches this group's base (if any).
  attach_sidecar <- function(group) {
    i <- which(sidecar_bases == group$base)
    if (length(i) >= 1) {
      group$metadata <- sidecars[i[1]]
      group$files <- c(group$files, sidecars[i[1]])
    } else {
      group$metadata <- NULL
    }
    group
  }

  exts <- tolower(tools::file_ext(data_files))
  bases <- tools::file_path_sans_ext(data_files)

  if (format %in% c("csv", "rds")) {
    allowed <- if (format == "rds") "rds" else c("csv", "txt", "tsv")
    idx <- which(exts %in% allowed)
    return(lapply(idx, function(i) {
      attach_sidecar(list(base = bases[i], files = data_files[i],
                          primary = data_files[i]))
    }))
  }

  primary_ext <- if (format == "imzML") "imzml" else "hdr"
  unique_bases <- unique(bases)
  lapply(unique_bases, function(b) {
    members <- data_files[bases == b]
    member_exts <- tolower(tools::file_ext(members))
    primary <- members[member_exts == primary_ext]
    attach_sidecar(list(
      base = b,
      files = members,
      primary = if (length(primary) >= 1) primary[1] else NA_character_
    ))
  })
}

# Restore added pixelData columns from a metadata sidecar CSV, matching rows
# by x/y/z/run (positional fallback when row counts agree).
apply_metadata_sidecar <- function(object, sidecar_path) {
  if (is.null(sidecar_path) || !file.exists(sidecar_path)) return(object)
  meta <- utils::read.csv(sidecar_path, check.names = FALSE,
                          stringsAsFactors = FALSE)
  if (!is.data.frame(meta) || nrow(meta) == 0) return(object)

  pd <- as.data.frame(Cardinal::pData(object))
  n <- nrow(pd)

  # Shared coordinate/run columns are identifiers; the rest are restored.
  id_cols <- intersect(names(meta), names(pd))
  added <- setdiff(names(meta), id_cols)
  if (length(added) == 0) return(object)

  order_idx <- NULL
  if (length(id_cols) > 0 && nrow(meta) >= n) {
    key_obj <- do.call(paste, c(lapply(id_cols, function(c) as.character(pd[[c]])),
                                sep = "\r"))
    key_meta <- do.call(paste, c(lapply(id_cols, function(c) as.character(meta[[c]])),
                                 sep = "\r"))
    m <- match(key_obj, key_meta)
    if (!anyNA(m)) order_idx <- m
  }
  if (is.null(order_idx)) {
    if (nrow(meta) != n) {
      stop(sprintf(paste("Metadata sidecar has %d rows but the dataset has %d",
                         "pixels and coordinates do not match; cannot reconcile."),
                   nrow(meta), n), call. = FALSE)
    }
    order_idx <- seq_len(n)
  }

  for (col in added) {
    Cardinal::pData(object)[[col]] <- meta[[col]][order_idx]
  }
  object
}

# Recover x/y from `px_<x>_<y>[_<z>]` pixel column names; NULL if any name
# does not match the convention.
parse_pixel_coord_names <- function(nms) {
  if (length(nms) == 0 || !all(grepl("^px_", nms))) return(NULL)
  parts <- strsplit(sub("^px_", "", nms), "_", fixed = TRUE)
  lens <- lengths(parts)
  if (length(unique(lens)) != 1 || lens[1] < 2) return(NULL)
  x <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 1)))
  y <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))
  if (anyNA(x) || anyNA(y)) return(NULL)
  data.frame(x = x, y = y)
}

# Detect the exporter's feature-by-pixel layout (mz column + px_* columns)
# and transpose it to one row per pixel; other layouts pass through.
normalize_csv_orientation <- function(df) {
  if (!is.data.frame(df) || ncol(df) < 2) return(df)
  nms <- names(df)
  mz_col <- which(tolower(nms) == "mz")[1]
  if (is.na(mz_col)) return(df)
  other <- setdiff(seq_along(nms), mz_col)
  # A numeric non-mz header means the frame is already pixel-by-feature.
  if (any(!is.na(suppressWarnings(as.numeric(nms[other]))))) return(df)
  mzvals <- suppressWarnings(as.numeric(df[[mz_col]]))
  if (anyNA(mzvals)) return(df)

  intensity <- t(as.matrix(df[, other, drop = FALSE]))  # pixels x features
  out <- as.data.frame(intensity, check.names = FALSE)
  names(out) <- as.character(mzvals)

  coords <- parse_pixel_coord_names(nms[other])
  if (!is.null(coords)) out <- cbind(coords, out)
  out
}

# One row per pixel; optional x/y/run columns (case-insensitive); remaining
# headers are numeric m/z values. 1-D layout when x/y absent.
parse_csv_msi <- function(df) {
  df <- normalize_csv_orientation(df)
  if (!is.data.frame(df) || ncol(df) < 2 || nrow(df) < 1) {
    stop("CSV must have at least one row and two columns.", call. = FALSE)
  }
  nms <- names(df)
  lower <- tolower(nms)
  x_col <- which(lower == "x")[1]
  y_col <- which(lower == "y")[1]
  run_col <- which(lower == "run")[1]
  meta_idx <- c(x_col, y_col, run_col)
  meta_idx <- meta_idx[!is.na(meta_idx)]

  feat_idx <- setdiff(seq_along(nms), meta_idx)
  mz <- suppressWarnings(as.numeric(nms[feat_idx]))
  ok <- !is.na(mz)
  if (sum(ok) < 1) {
    stop("No numeric m/z columns detected. Intensity column headers must be m/z values.",
         call. = FALSE)
  }
  feat_idx <- feat_idx[ok]
  mz <- mz[ok]
  ord <- order(mz)
  mz <- mz[ord]
  feat_idx <- feat_idx[ord]

  intensity <- t(as.matrix(df[, feat_idx, drop = FALSE]))
  storage.mode(intensity) <- "double"
  dimnames(intensity) <- NULL

  n_pix <- nrow(df)
  if (!is.na(x_col) && !is.na(y_col)) {
    coord <- data.frame(x = as.numeric(df[[x_col]]), y = as.numeric(df[[y_col]]))
  } else {
    coord <- data.frame(x = seq_len(n_pix), y = rep(1, n_pix))
  }
  run <- if (!is.na(run_col)) {
    factor(as.character(df[[run_col]]))
  } else {
    factor(rep("csv", n_pix))
  }

  list(mz = mz, coord = coord, run = run, intensity = intensity,
       n_pixels = n_pix, n_features = length(mz))
}

# Rename a dataset record; errors on unknown id or empty name.
rename_dataset_record <- function(store, id, new_name) {
  if (is.null(id) || !nzchar(id) || is.null(store[[id]])) {
    stop("Unknown dataset id.", call. = FALSE)
  }
  nm <- trimws(new_name %||% "")
  if (!nzchar(nm)) {
    stop("New name must not be empty.", call. = FALSE)
  }
  store[[id]]$name <- nm
  store
}

build_msi_from_csv <- function(df) {
  parts <- parse_csv_msi(df)
  fd <- Cardinal::MassDataFrame(mz = parts$mz)
  pd <- Cardinal::PositionDataFrame(coord = parts$coord, run = parts$run)
  Cardinal::MSImagingExperiment(spectraData = parts$intensity,
                                featureData = fd, pixelData = pd)
}

# CSV/TSV/TXT to MSImagingExperiment; sep defaults to tab for .tsv.
read_csv_msi <- function(path, sep = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (is.null(sep)) sep <- if (ext == "tsv") "\t" else ","
  df <- utils::read.csv(path, sep = sep, check.names = FALSE,
                        stringsAsFactors = FALSE)
  build_msi_from_csv(df)
}

# .rds to MSImagingExperiment; MSImagingArrays (per-pixel m/z axes) are
# coerced for consistency.
read_rds_msi <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  obj <- tryCatch(readRDS(path), error = function(e) {
    stop(sprintf("Could not read RDS file: %s", conditionMessage(e)),
         call. = FALSE)
  })
  if (!is(obj, "MSImagingExperiment") && !is(obj, "MSImagingArrays")) {
    stop(sprintf(paste("RDS file does not contain a Cardinal MSI object",
                       "(found '%s'). Expected MSImagingExperiment."),
                 class(obj)[1]), call. = FALSE)
  }
  if (!is(obj, "MSImagingExperiment")) {
    obj <- as(obj, "MSImagingExperiment")
  }
  obj
}

# Read an imzML with memory mode, mass-range, and resolution options.
read_imzml <- function(path, memory = "memory", mass_range = NULL,
                       resolution = NULL, units = c("ppm", "mz")) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  units <- match.arg(units)
  memory <- match.arg(memory, c("memory", "shared"))
  params <- list(file = path, memory = memory)
  if (!is.null(mass_range) && !all(is.na(mass_range))) {
    if (length(mass_range) != 2 || anyNA(mass_range)) {
      stop("mass_range must be a length-2 numeric vector c(min, max).",
           call. = FALSE)
    }
    if (any(mass_range < 0)) {
      stop("mass_range values must be non-negative.", call. = FALSE)
    }
    if (mass_range[1] >= mass_range[2]) {
      stop("mass_range minimum must be less than the maximum.", call. = FALSE)
    }
    params$mass.range <- mass_range
  }
  if (!is.null(resolution) && !is.na(resolution)) {
    if (resolution <= 0) {
      stop("resolution must be a positive number.", call. = FALSE)
    }
    params$resolution <- resolution
    # units: "ppm" (spacing grows with m/z) or "mz" (fixed Da spacing).
    params$units <- units
  }
  obj <- do.call(Cardinal::readMSIData, params)

  # Processed imzML may return MSImagingArrays (per-pixel m/z axes); the app
  # assumes a common feature axis, so coerce first.
  if (!is(obj, "MSImagingExperiment")) {
    obj <- as(obj, "MSImagingExperiment")
  }

  # readMSIData honours mass.range only when resampling a processed imzML;
  # continuous imzML ignores it, so apply an explicit feature subset.
  if (!is.null(params$mass.range)) {
    mzv <- Cardinal::mz(obj)
    keep <- !is.na(mzv) & mzv >= mass_range[1] & mzv <= mass_range[2]
    if (!any(keep)) {
      stop(sprintf("No m/z values fall within the mass range %g-%g.",
                   mass_range[1], mass_range[2]), call. = FALSE)
    }
    if (!all(keep)) obj <- obj[keep, ]
  }
  obj
}

# Normalized imzML options list: memory, mass_range, resolution, units.
build_imzml_opts <- function(memory = "memory", mass_min = NA, mass_max = NA,
                             resolution = NA, units = "ppm") {
  mass_range <- NULL
  if (!is.na(mass_min) && !is.na(mass_max)) {
    mass_range <- c(mass_min, mass_max)
  }
  res <- if (!is.na(resolution)) resolution else NULL
  units <- match.arg(units, c("ppm", "mz"))
  list(memory = memory, mass_range = mass_range, resolution = res,
       units = units)
}

# Importable file names (not paths) under `dir`: imzml/ibd, hdr/img/t2m,
# csv/txt/tsv, or rds by format.
list_import_files <- function(dir, format, recursive = FALSE) {
  if (!dir.exists(dir)) {
    stop(sprintf("Folder not found: %s", dir), call. = FALSE)
  }
  allowed <- switch(format,
                    imzML = c("imzml", "ibd"),
                    analyze = c("hdr", "img", "t2m"),
                    csv = c("csv", "txt", "tsv"),
                    rds = "rds",
                    stop(sprintf("Unknown format '%s'.", format), call. = FALSE))
  all_files <- list.files(dir, recursive = recursive, full.names = FALSE)
  exts <- tolower(tools::file_ext(all_files))
  keep <- exts %in% allowed
  # Also surface metadata sidecars (rds is self-contained, skips them).
  if (!identical(format, "rds")) {
    keep <- keep | grepl(import_metadata_pattern, all_files, ignore.case = TRUE)
  }
  all_files[keep]
}

dataImportUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
        h3("Data Import"),
        p("Import mass spectrometry imaging data in imzML, Analyze 7.5, CSV, or Rds format")
      )
    ),
    
    fluidRow(
      column(6,
        bslib::card(
          bslib::card_header("File Selection"),
          bslib::card_body(
          
          selectInput(ns("format"), "File Format",
                      choices = c("imzML" = "imzML",
                                  "Analyze 7.5" = "analyze",
                                  "CSV/Tabular" = "csv",
                                  "R data (.rds)" = "rds"),
                      selected = "imzML"),
          
          uiOutput(ns("file_input_label")),
          fileInput(ns("file"), "Choose File(s)",
                    accept = NULL,
                    multiple = TRUE),
          
          actionButton(ns("import_btn"), "Import Data", 
                      class = "btn-primary"),

          tags$hr(),
          tags$label("Bulk import from a folder"),
          helpText(paste("Import every matching file in a server-side folder in one",
                         "action. Files are auto-grouped into datasets by base name",
                         "(imzML/ibd, hdr/img/t2m).")),
          textInput(ns("folder_path"), "Folder path",
                    placeholder = "e.g. C:/data/msi_batch"),
          checkboxInput(ns("folder_recursive"), "Include subfolders", value = FALSE),
          actionButton(ns("import_folder_btn"), "Import Folder",
                      icon = icon("folder-open"), class = "btn-outline-primary")
          )
        )
      ),
      
      column(6,
        bslib::card(
          bslib::card_header("Advanced Settings"),
          bslib::card_body(
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'imzML'", ns("format")),
            selectInput(ns("memory_mode"), "Memory mode",
                       choices = c("In-memory (faster)" = "memory",
                                   "File-backed (lower memory)" = "shared"),
                       selected = "memory"),
            helpText("File-backed mode keeps spectra on disk, using less RAM but slower access. Recommended for datasets larger than available memory."),
            
            numericInput(ns("mass_min"), "Minimum m/z",
                        value = NA, min = 0, step = 0.1),
            numericInput(ns("mass_max"), "Maximum m/z",
                        value = NA, min = 0, step = 0.1),
            numericInput(ns("resolution"), "Resolution",
                        value = NA, min = 0, step = 0.01),
            selectInput(ns("resolution_units"), "Resolution units",
                        choices = c("ppm" = "ppm",
                                    "mz" = "mz"),
                        selected = "ppm"),
            helpText(paste("Resolution defines the common m/z axis spacing when a",
                           "processed imzML is resampled. Units are only applied",
                           "when a resolution is set."))
          ),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'analyze'", ns("format")),
            helpText("Advanced settings are only available for imzML format.")
          ),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'csv'", ns("format")),
            helpText(strong("Expected CSV/tabular layout:")),
            tags$ul(
              tags$li("One row per pixel/spectrum."),
              tags$li("Optional 'x', 'y', and 'run' columns (case-insensitive)."),
              tags$li("Every other column header must be a numeric m/z value."),
              tags$li("Without 'x'/'y', pixels are placed in a 1-D layout.")
            ),
            helpText("Multiple tabular files import as separate datasets.")
          ),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'rds'", ns("format")),
            helpText(strong("R data (.rds) import:")),
            tags$ul(
              tags$li("The file must contain a Cardinal MSI object saved with the Export tab (or saveRDS)."),
              tags$li("Each .rds file imports as its own dataset.")
            ),
            helpText("Advanced settings do not apply to .rds files.")
          )
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        uiOutput(ns("import_status")),
        uiOutput(ns("data_summary"))
      )
    ),
    
    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header("Imported Datasets"),
          bslib::card_body(
          p("Select a row to make that dataset active. Use the table headers to sort or filter."),
          DT::dataTableOutput(ns("dataset_table")),
          uiOutput(ns("active_dataset_label")),
          tags$hr(),
          tags$strong("Selected dataset actions"),
          fluidRow(
            column(8,
              textInput(ns("rename_value"), "Rename to",
                        placeholder = "New name for selected dataset",
                        width = "100%")
            )
          ),
          div(class = "d-flex gap-2",
            actionButton(ns("rename_dataset"), "Rename",
                         class = "btn-primary"),
            actionButton(ns("delete_dataset"), "Delete",
                         class = "btn-danger")
          )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header("Data Preview"),
          bslib::card_body(
            uiOutput(ns("preview_header")),
            fluidRow(
              column(7,
                shinycssloaders::withSpinner(
                  plotOutput(ns("preview_plot"), height = "300px")
                )
              ),
              column(5,
                strong("First pixels"),
                DT::dataTableOutput(ns("preview_table"))
              )
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        bslib::card(
          bslib::card_header("Import History"),
          bslib::card_body(
            p("Chronological log of all import attempts this session."),
            DT::dataTableOutput(ns("history_table"))
          )
        )
      )
    )
  )
}

dataImportServer <- function(id, datasets, active_id, perf_log = NULL) {
  moduleServer(id, function(input, output, session) {
    
    library(Cardinal)
    library(DT)

    # File-backed datasets keep a live handle to their backing files; remove
    # per-import temp dirs only when the session ends.
    import_temp_dirs <- character(0)
    session$onSessionEnded(function() {
      for (d in import_temp_dirs) unlink(d, recursive = TRUE)
    })

    import_status <- reactiveVal(list(
      status = "idle",
      message = "Ready to import data",
      type = "info"
    ))
    
    dataset_counter <- reactiveVal(0)

    # Chronological log of import attempts (most recent first)
    import_history <- reactiveVal(list())

    log_history <- function(name, format, status, message,
                            pixels = NA, features = NA) {
      entry <- list(
        time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        name = name,
        format = format,
        status = status,
        pixels = pixels,
        features = features,
        message = message
      )
      import_history(c(list(entry), import_history()))
    }
    
    output$file_input_label <- renderUI({
      label_text <- switch(input$format,
                           "imzML" = "Choose imzML file (and .ibd if needed)",
                           "analyze" = "Choose Analyze 7.5 files (.hdr, .img, .t2m)",
                           "csv" = "Choose CSV/Tabular file",
                           "rds" = "Choose R data (.rds) file")
      strong(label_text)
    })
    
    observeEvent(input$file, {
      req(input$file)
      
      # Metadata sidecars are not datasets; exclude from pairing validation.
      files <- input$file
      all_names <- files$name
      is_sidecar <- grepl(import_metadata_pattern, all_names, ignore.case = TRUE)
      file_names <- all_names[!is_sidecar]
      file_exts <- tools::file_ext(file_names)
      file_bases <- tools::file_path_sans_ext(file_names)
      
      if (input$format == "imzML") {
        unique_bases <- unique(file_bases)
        valid_pairs <- TRUE
        missing_pairs <- c()
        
        for (base in unique_bases) {
          has_imzml <- any(grepl(paste0("^", base, "\\.imzML$"), file_names, ignore.case = TRUE))
          has_ibd <- any(grepl(paste0("^", base, "\\.ibd$"), file_names, ignore.case = TRUE))
          
          if (!has_imzml || !has_ibd) {
            valid_pairs <- FALSE
            missing_pairs <- c(missing_pairs, base)
          }
        }
        
        if (!valid_pairs) {
          import_status(list(
            status = "error",
            message = paste("Missing paired files for:", paste(missing_pairs, collapse = ", "), 
                          "\nPlease select both .imzML and .ibd files for each dataset."),
            type = "danger"
          ))
        } else {
          import_status(list(
            status = "ready",
            message = sprintf("Valid file pairs detected: %d dataset(s). Click 'Import Data' to proceed.", 
                            length(unique_bases)),
            type = "info"
          ))
        }
      } else if (input$format == "analyze") {
        unique_bases <- unique(file_bases)
        valid_pairs <- TRUE
        missing_pairs <- c()
        
        for (base in unique_bases) {
          has_hdr <- any(grepl(paste0("^", base, "\\.hdr$"), file_names, ignore.case = TRUE))
          has_img <- any(grepl(paste0("^", base, "\\.img$"), file_names, ignore.case = TRUE))
          has_t2m <- any(grepl(paste0("^", base, "\\.t2m$"), file_names, ignore.case = TRUE))
          
          if (!has_hdr || !has_img || !has_t2m) {
            valid_pairs <- FALSE
            missing_pairs <- c(missing_pairs, base)
          }
        }
        
        if (!valid_pairs) {
          import_status(list(
            status = "error",
            message = paste("Missing paired files for:", paste(missing_pairs, collapse = ", "),
                          "\nPlease select .hdr, .img, and .t2m files for each Analyze 7.5 dataset."),
            type = "danger"
          ))
        } else {
          import_status(list(
            status = "ready",
            message = sprintf("Valid file pairs detected: %d dataset(s). Click 'Import Data' to proceed.",
                            length(unique_bases)),
            type = "info"
          ))
        }
      } else if (input$format == "csv") {
        if (!all(file_exts %in% c("csv", "txt", "tsv"))) {
          import_status(list(
            status = "error",
            message = "Invalid file extension for CSV format. Expected .csv, .txt, or .tsv",
            type = "danger"
          ))
        } else {
          import_status(list(
            status = "ready",
            message = sprintf("%d file(s) selected. Click 'Import Data' to proceed.", length(file_names)),
            type = "info"
          ))
        }
      } else if (input$format == "rds") {
        if (!all(tolower(file_exts) == "rds")) {
          import_status(list(
            status = "error",
            message = "Invalid file extension for R data format. Expected .rds",
            type = "danger"
          ))
        } else {
          import_status(list(
            status = "ready",
            message = sprintf("%d file(s) selected. Click 'Import Data' to proceed.", length(file_names)),
            type = "info"
          ))
        }
      }
    })
    
    # Normalized imzML options from the advanced-settings controls.
    imzml_import_opts <- reactive({
      build_imzml_opts(memory = input$memory_mode,
                       mass_min = input$mass_min,
                       mass_max = input$mass_max,
                       resolution = input$resolution,
                       units = input$resolution_units %||% "ppm")
    })

    # Read one dataset group; imzml_opts ignored for other formats.
    read_group <- function(group, temp_dir, format, imzml_opts = NULL) {
      data <- if (format == "csv") {
        read_csv_msi(file.path(temp_dir, group$primary))
      } else if (format == "rds") {
        read_rds_msi(file.path(temp_dir, group$primary))
      } else if (is.na(group$primary)) {
        stop(sprintf("No primary file found for '%s'.", group$base), call. = FALSE)
      } else if (format == "imzML") {
        if (is.null(imzml_opts)) imzml_opts <- build_imzml_opts()
        read_imzml(file.path(temp_dir, group$primary),
                   memory = imzml_opts$memory,
                   mass_range = imzml_opts$mass_range,
                   resolution = imzml_opts$resolution,
                   units = imzml_opts$units %||% "ppm")
      } else {
        Cardinal::readMSIData(file.path(temp_dir, group$primary))
      }

      # Restore added pixelData columns from a matching sidecar, if any.
      if (!is.null(group$metadata)) {
        data <- apply_metadata_sidecar(data, file.path(temp_dir, group$metadata))
      }
      data
    }

    # Import each group from `source_dir`; shared by upload/folder observers.
    import_groups <- function(groups, source_dir, fmt) {
      if (length(groups) == 0) {
        import_status(list(status = "error",
                           message = "No valid files detected for the selected format.",
                           type = "danger"))
        return()
      }

      n_ok <- 0
      n_fail <- 0
      last_id <- NULL
      import_t0 <- proc.time()[["elapsed"]]

      withProgress(message = "Importing datasets", value = 0, {
        n <- length(groups)
        for (gi in seq_len(n)) {
          group <- groups[[gi]]
          # Show only completed work; jumping ahead hit 100% while reading.
          setProgress(value = (gi - 1) / n,
                      detail = sprintf("Reading (%d/%d) %s", gi, n, group$base))

          tryCatch({
            data <- read_group(group, source_dir, fmt,
                               imzml_opts = imzml_import_opts())
            summary <- summarize_msi(data)

            # Id unique across the store, incl. save-flow datasets.
            new_id <- next_dataset_id(names(datasets()))
            dataset_counter(dataset_counter() + 1)

            record <- list(
              id = new_id,
              name = group$base,
              format = fmt,
              files = paste(group$files, collapse = ", "),
              extensions = paste(unique(tolower(tools::file_ext(group$files))),
                                 collapse = ", "),
              imported_at = Sys.time(),
              summary = summary,
              data = data
            )

            current <- datasets()
            current[[new_id]] <- record
            datasets(current)
            last_id <- new_id
            n_ok <- n_ok + 1

            log_history(group$base, fmt, "success",
                        sprintf("%d pixels, %d features",
                                summary$n_pixels, summary$n_features),
                        pixels = summary$n_pixels, features = summary$n_features)
          }, error = function(e) {
            n_fail <<- n_fail + 1
            log_history(group$base, fmt, "failed", conditionMessage(e))
          })

          setProgress(value = gi / n,
                      detail = sprintf("(%d/%d) done", gi, n))
        }
      })

      if (!is.null(last_id)) {
        active_id(last_id)
      }

      if (!is.null(perf_log)) {
        perf_log(perf_log_append(perf_log(), "import",
                                 sprintf("import %s (%d ok, %d failed)",
                                         fmt, n_ok, n_fail),
                                 proc.time()[["elapsed"]] - import_t0))
      }

      if (n_ok > 0 && n_fail == 0) {
        import_status(list(status = "success",
                           message = sprintf("Successfully imported %d dataset(s).", n_ok),
                           type = "success"))
      } else if (n_ok > 0 && n_fail > 0) {
        import_status(list(status = "warning",
                           message = sprintf("Imported %d dataset(s); %d failed. See history below.",
                                             n_ok, n_fail),
                           type = "warning"))
      } else {
        import_status(list(status = "error",
                           message = sprintf("All %d import(s) failed. See history below.", n_fail),
                           type = "danger"))
      }
    }

    observeEvent(input$import_btn, {
      req(input$file)

      files <- input$file
      file_names <- files$name
      fmt <- input$format

      # Move uploads into one temp dir so Cardinal can resolve sidecars by
      # name; fall back to copy on failure (e.g. across filesystems).
      temp_dir <- tempfile(pattern = "cardinal_nest_")
      dir.create(temp_dir)
      for (i in seq_len(nrow(files))) {
        dest <- file.path(temp_dir, files$name[i])
        if (!file.rename(files$datapath[i], dest)) {
          file.copy(files$datapath[i], dest)
        }
      }
      # Keep backing files alive for the session (file-backed datasets).
      import_temp_dirs <<- c(import_temp_dirs, temp_dir)

      groups <- import_group_files(file_names, fmt)
      import_groups(groups, temp_dir, fmt)

      # Datapaths are stale after the move; reset so the same file re-imports.
      shinyjs::reset("file")
    })

    # Bulk import from a server-side folder (no copy needed).
    observeEvent(input$import_folder_btn, {
      folder <- trimws(input$folder_path %||% "")
      fmt <- input$format
      if (!nzchar(folder)) {
        import_status(list(status = "error",
                           message = "Enter a folder path to import from.",
                           type = "danger"))
        return()
      }

      file_names <- tryCatch(
        list_import_files(folder, fmt, recursive = isTRUE(input$folder_recursive)),
        error = function(e) {
          import_status(list(status = "error", message = conditionMessage(e),
                             type = "danger"))
          NULL
        }
      )
      if (is.null(file_names)) return()

      groups <- import_group_files(file_names, fmt)
      import_groups(groups, folder, fmt)
    })
    
    output$import_status <- render_status(import_status)
    
    output$data_summary <- renderUI({
      ds <- datasets()
      id <- active_id()
      if (is.null(id) || is.null(ds[[id]])) {
        return(NULL)
      }
      
      record <- ds[[id]]
      summary <- record$summary
      
      tagList(
        h4(sprintf("Data Summary: %s", record$name)),
        fluidRow(
          column(4,
            strong("Pixels:"), summary$n_pixels, br(),
            strong("Features:"), summary$n_features, br(),
            strong("Runs:"), summary$n_runs
          ),
          column(4,
            strong("m/z Range:"), 
            sprintf("%.2f - %.2f", summary$mz_range[1], summary$mz_range[2]), br(),
            strong("Centroided:"), summary$centroided, br(),
            strong("Class:"), summary$class
          ),
          column(4,
            strong("Format:"), record$format
          )
        )
      )
    })
    
    output$preview_header <- renderUI({
      ds <- datasets()
      id <- active_id()
      if (is.null(id) || is.null(ds[[id]])) {
        return(helpText("Import or select a dataset to preview it."))
      }
      helpText(strong("Previewing: "), ds[[id]]$name,
               " (mean spectrum across all pixels)")
    })

    output$preview_plot <- renderPlot({
      ds <- datasets()
      id <- active_id()
      if (is.null(id) || is.null(ds[[id]])) return(NULL)
      data <- ds[[id]]$data

      mzv <- Cardinal::mz(data)
      mean_int <- rowMeans(as.matrix(Cardinal::spectra(data)), na.rm = TRUE)
      plot(mzv, mean_int, type = "l",
           xlab = "m/z", ylab = "Mean intensity",
           main = "Mean spectrum")
    })

    output$preview_table <- DT::renderDataTable({
      ds <- datasets()
      id <- active_id()
      if (is.null(id) || is.null(ds[[id]])) {
        return(DT::datatable(data.frame(Message = "No active dataset"),
                             options = list(dom = "t"), rownames = FALSE))
      }
      cd <- as.data.frame(Cardinal::coord(ds[[id]]$data))
      DT::datatable(utils::head(cd, 10),
                    options = list(dom = "t", scrollX = TRUE),
                    rownames = FALSE)
    })

    output$history_table <- DT::renderDataTable({
      hist <- import_history()
      if (length(hist) == 0) {
        return(DT::datatable(data.frame(Message = "No imports yet this session"),
                             options = list(dom = "t"), rownames = FALSE))
      }
      df <- data.frame(
        Time = vapply(hist, function(h) h$time, character(1)),
        Name = vapply(hist, function(h) h$name, character(1)),
        Format = vapply(hist, function(h) h$format, character(1)),
        Status = vapply(hist, function(h) h$status, character(1)),
        Pixels = vapply(hist, function(h) as.character(h$pixels), character(1)),
        Features = vapply(hist, function(h) as.character(h$features), character(1)),
        Details = vapply(hist, function(h) h$message, character(1)),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, options = list(pageLength = 5, scrollX = TRUE),
                    rownames = FALSE)
    })

    # Stable ordering: table row index maps to dataset id.
    dataset_overview <- reactive({
      ds <- datasets()
      if (length(ds) == 0) return(NULL)
      
      data.frame(
        Name = vapply(ds, function(r) r$name, character(1)),
        Format = vapply(ds, function(r) r$format, character(1)),
        Extensions = vapply(ds, function(r) r$extensions, character(1)),
        Files = vapply(ds, function(r) r$files, character(1)),
        Pixels = vapply(ds, function(r) r$summary$n_pixels, numeric(1)),
        Features = vapply(ds, function(r) r$summary$n_features, numeric(1)),
        Imported = vapply(ds, function(r) format(r$imported_at, "%Y-%m-%d %H:%M:%S"),
                          character(1)),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })
    
    output$dataset_table <- DT::renderDataTable({
      overview <- dataset_overview()
      if (is.null(overview)) {
        return(DT::datatable(
          data.frame(Message = "No datasets imported yet"),
          options = list(dom = "t"),
          rownames = FALSE
        ))
      }
      
      ids <- names(datasets())
      selected_row <- which(ids == active_id())
      
      DT::datatable(
        overview,
        selection = list(mode = "single", selected = selected_row),
        options = list(pageLength = 5, scrollX = TRUE),
        filter = "top",
        rownames = FALSE
      )
    })
    
    observeEvent(input$dataset_table_rows_selected, {
      selected <- input$dataset_table_rows_selected
      ids <- names(datasets())
      if (length(selected) == 1 && selected <= length(ids)) {
        active_id(ids[selected])
      }
    })
    
    observeEvent(input$rename_dataset, {
      selected <- input$dataset_table_rows_selected
      ids <- names(datasets())

      if (length(selected) != 1 || selected > length(ids)) {
        showNotification("Select a dataset row to rename", type = "warning")
        return()
      }

      new_name <- trimws(input$rename_value %||% "")
      if (!nzchar(new_name)) {
        showNotification("Enter a new name before renaming", type = "warning")
        return()
      }

      rid <- ids[selected]
      updated <- tryCatch(
        rename_dataset_record(datasets(), rid, new_name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error")
          NULL
        }
      )
      if (is.null(updated)) return()

      datasets(updated)
      updateTextInput(session, "rename_value", value = "")
      showNotification(sprintf("Dataset renamed to '%s'", new_name),
                       type = "message")
    })

    observeEvent(input$delete_dataset, {
      selected <- input$dataset_table_rows_selected
      ids <- names(datasets())
      
      if (length(selected) != 1 || selected > length(ids)) {
        showNotification("Select a dataset row to delete", type = "warning")
        return()
      }
      
      del_id <- ids[selected]
      current <- datasets()
      current[[del_id]] <- NULL
      datasets(current)
      
      # Reassign active dataset if the deleted one was active.
      remaining <- names(current)
      if (identical(active_id(), del_id)) {
        active_id(if (length(remaining) > 0) remaining[length(remaining)] else NULL)
      }
      
      showNotification("Dataset deleted", type = "message")
    })
    
    output$active_dataset_label <- renderUI({
      ds <- datasets()
      id <- active_id()
      if (is.null(id) || is.null(ds[[id]])) {
        return(helpText("No active dataset."))
      }
      helpText(strong("Active dataset: "), ds[[id]]$name)
    })
  })
}
