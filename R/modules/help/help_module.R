# Cardinal NEST - help module: log viewer and vignette browser

library(shiny)

# Remove a leading YAML frontmatter block, if present.
strip_yaml_frontmatter <- function(lines) {
  if (length(lines) == 0) return(lines)
  if (!grepl("^---\\s*$", lines[1])) return(lines)
  end <- grep("^---\\s*$", lines)[2]
  if (is.na(end)) return(lines)
  lines[-seq_len(end)]
}

# "01-installation-and-first.Rmd" -> "Installation And First".
make_vignette_title <- function(path) {
  nm <- tools::file_path_sans_ext(basename(path))
  nm <- sub("^[0-9]+[-_]", "", nm)
  nm <- gsub("[-_]", " ", nm)
  tools::toTitleCase(nm)
}

# Convert Rmd code chunks to plain markdown fences; chunks with
# include=FALSE or echo=FALSE are dropped entirely.
prepare_vignette_md <- function(lines) {
  lines <- strip_yaml_frontmatter(lines)
  out <- character(0)
  i <- 1
  while (i <= length(lines)) {
    if (grepl("^```\\{", lines[i])) {
      m <- regexec("^```\\{([^}]*)\\}", lines[i])
      cfg <- regmatches(lines[i], m)[[1]][2]
      drop <- grepl("\\b(include|echo)\\s*=\\s*FALSE", cfg, ignore.case = TRUE)
      j <- i + 1
      while (j <= length(lines) && !grepl("^```\\s*$", lines[j])) j <- j + 1
      if (!drop) {
        lang <- "r"
        has_engine <- grepl("engine\\s*=\\s*['\"]", cfg)
        if (has_engine) {
          lang <- sub(".*engine\\s*=\\s*['\"]([^'\"]+)['\"].*", "\\1", cfg)
        } else {
          first <- sub("\\s.*", "", cfg)
          if (nzchar(first) && !grepl("=", first)) lang <- first
        }
        if (i < j - 1) {
          out <- c(out, sprintf("```%s", lang), lines[(i + 1):(j - 1)], "```")
        } else {
          out <- c(out, sprintf("```%s", lang), "```")
        }
      }
      i <- j + 1
    } else {
      out <- c(out, lines[i])
      i <- i + 1
    }
  }
  out
}

# Markdown to HTML via commonmark; falls back to escaped preformatted text.
render_vignette_html <- function(md_text) {
  if (requireNamespace("commonmark", quietly = TRUE)) {
    commonmark::markdown_html(md_text, extensions = TRUE)
  } else {
    as.character(tags$pre(style = "white-space: pre-wrap;", md_text))
  }
}

helpUI <- function(id) {
  ns <- NS(id)

  tagList(
    h3("Help & Logs"),
    bslib::navset_tab(
      bslib::nav_panel("Log",
        p(class = "text-muted small",
          "Console messages, warnings, and progress output from the current ",
          "session are captured here in real time."),
        verbatimTextOutput(ns("log"), placeholder = TRUE)
      ),
      bslib::nav_panel("Vignette",
        p(class = "text-muted small",
          "Select a vignette to view the user guide inside the app."),
        selectInput(ns("vignette_select"), "Vignette", choices = NULL),
        htmlOutput(ns("vignette_view"))
      )
    )
  )
}

helpServer <- function(id, log_path,
                       vignette_dir = file.path(getwd(), "docs", "tutorials")) {
  moduleServer(id, function(input, output, session) {
    log_lines <- reactivePoll(
      intervalMillis = 1000,
      session = session,
      checkFunc = function() {
        if (file.exists(log_path)) file.info(log_path)$mtime else 0
      },
      valueFunc = function() {
        if (!file.exists(log_path)) return(character(0))
        readLines(log_path, warn = FALSE)
      }
    )

    output$log <- renderText({
      lines <- log_lines()
      if (length(lines) > 500) lines <- utils::tail(lines, 500)
      paste(lines, collapse = "\n")
    })

    observe({
      if (!dir.exists(vignette_dir)) {
        updateSelectInput(session, "vignette_select",
                          choices = c("No vignettes found" = ""))
        return()
      }
      files <- sort(list.files(vignette_dir, pattern = "\\.Rmd$",
                               full.names = TRUE, ignore.case = TRUE))
      if (length(files) == 0) {
        updateSelectInput(session, "vignette_select",
                          choices = c("No vignettes found" = ""))
        return()
      }
      choices <- stats::setNames(basename(files), vapply(files, make_vignette_title, character(1)))
      selected <- if (length(input$vignette_select) &&
                       input$vignette_select %in% choices) {
        input$vignette_select
      } else {
        default <- grep("installation.*first.*launch", basename(files),
                        value = TRUE, ignore.case = TRUE)[1]
        if (is.na(default)) default <- choices[1]
        default
      }
      updateSelectInput(session, "vignette_select", choices = choices, selected = selected)
    })

    vignette_html <- reactive({
      req(input$vignette_select)
      req(nzchar(input$vignette_select))
      path <- file.path(vignette_dir, input$vignette_select)
      req(file.exists(path))
      lines <- readLines(path, warn = FALSE)
      lines <- prepare_vignette_md(lines)
      body <- paste(lines, collapse = "\n")
      render_vignette_html(body)
    })

    output$vignette_view <- renderUI({
      html <- vignette_html()
      tagList(
        tags$style(HTML(".vignette-content, .vignette-content pre, .vignette-content code { box-sizing: border-box; } .vignette-content { width: 100%; max-width: 100%; } .vignette-content pre, .vignette-content code { visibility: visible !important; } .vignette-content pre { width: 100%; max-width: 100%; display: block !important; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 0.25rem; padding: 0.75rem; color: #212529; font-family: SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; white-space: pre-wrap !important; word-wrap: break-word !important; overflow-wrap: anywhere !important; word-break: break-all !important; margin-bottom: 1rem; } .vignette-content code { background-color: rgba(0,0,0,0.08); color: #212529; padding: 0.15em 0.3em; border-radius: 0.2rem; font-family: SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; } .vignette-content pre > code { background-color: transparent; padding: 0; color: inherit; white-space: inherit !important; word-break: inherit !important; }")),
        div(class = "card p-3 vignette-content",
            style = "max-height: 210vh; overflow-y: auto;",
            HTML(html))
      )
    })
  })
}
