# Cardinal NEST - package entry point

# Tutorials ship in inst/tutorials; fall back to the source tree when the
# app runs from a cloned folder via app.R.
tutorials_dir <- function() {
  installed <- system.file("tutorials", package = "CardinalNEST")
  if (nzchar(installed)) return(installed)
  file.path(getwd(), "inst", "tutorials")
}

cardinal_nest_app <- function(log_file = file.path(getwd(), "logs", "app.log"),
                              options = list()) {
  log_file <- normalizePath(log_file, winslash = "/", mustWork = FALSE)

  # Capture console output/messages to a log file for the Help tab.
  start_logging <- function() {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    con <- file(log_file, open = "w")
    sink(con, type = "output", split = TRUE)
    sink(con, type = "message")
    onStop(function() {
      sink(type = "message")
      sink(type = "output")
      if (isOpen(con)) close(con)
    })
  }

  shinyApp(
    ui = main_ui(),
    server = function(input, output, session) {
      main_server(input, output, session,
                  log_path = log_file, vignette_dir = tutorials_dir())
    },
    onStart = function() {
      start_logging()
      # Default 30 MB upload limit is too small for imzML/ibd datasets.
      old_opts <- base::options(shiny.maxRequestSize = 10 * 1024^3)
      onStop(function() base::options(old_opts))
    },
    options = options
  )
}

#' Launch the Cardinal NEST Shiny application
#'
#' @param host Address to listen on. Use "0.0.0.0" to accept connections from
#'   other machines (for example a compute node reached through an SSH tunnel).
#' @param port Port to listen on.
#' @param launch.browser Open the app in a web browser.
#' @param log_file Path of the session log shown in the Help tab.
#' @param ... Further arguments passed to [shiny::runApp()].
#' @export
runCardinalNEST <- function(host = "127.0.0.1", port = 8080,
                            launch.browser = interactive(),
                            log_file = file.path(getwd(), "logs", "app.log"),
                            ...) {
  runApp(cardinal_nest_app(log_file = log_file),
         host = host, port = port, launch.browser = launch.browser, ...)
}
