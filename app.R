# Cardinal NEST - application entry point

library(shiny)
library(Cardinal)
library(bslib)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)

source("R/utils/helpers.R")
source("R/utils/palette.R")
source("R/modules/data_management/data_import_module.R")
source("R/modules/data_management/metadata_module.R")
source("R/modules/preprocessing/preprocessing_module.R")
source("R/modules/visualization/visualization_module.R")
source("R/modules/analysis/analysis_module.R")
source("R/modules/spatial/spatial_module.R")
source("R/modules/export/export_module.R")
source("R/modules/workflow/workflow_module.R")
source("R/modules/performance/performance_module.R")
source("R/modules/help/help_module.R")
source("R/ui/main_ui.R")
source("R/server/main_server.R")

# Capture console output/messages to a log file for the Help tab.
app_log_path <<- file.path(getwd(), "logs", "app.log")
dir.create(dirname(app_log_path), showWarnings = FALSE, recursive = TRUE)

app_log_con <<- file(app_log_path, open = "w")
sink(app_log_con, type = "output", split = TRUE)
sink(app_log_con, type = "message")

onStop(function() {
  sink(type = "message")
  sink(type = "output")
  if (isOpen(app_log_con)) close(app_log_con)
})

# Default 30 MB upload limit is too small for imzML/ibd datasets.
options(shiny.maxRequestSize = 10 * 1024^3)

shinyApp(
  ui = main_ui(),
  server = main_server,
  options = list(
    port = 8080,
    host = "127.0.0.1"
  )
)
