# Cardinal NEST - main server

library(shiny)

main_server <- function(input, output, session) {

  # Dataset store: named list of records (metadata + MSImagingExperiment).
  datasets <- reactiveVal(list())
  active_id <- reactiveVal(NULL)

  perf_log <- reactiveVal(NULL)
  perf_state <- reactiveValues(cache_enabled = TRUE)
  analysis_cache <- new_result_cache()

  mse_data <- active_dataset(datasets, active_id)

  dataImportServer("import", datasets, active_id, perf_log = perf_log)
  metadataServer("explore", datasets, active_id)

  # Returns the current pipeline as reusable specs for the Workflow tab.
  preproc <- preprocessingServer("preprocess", datasets, active_id,
                                 perf_log = perf_log)

  visualizationServer("visualize_content", datasets, active_id)

  analysisServer("analyze", datasets, active_id,
                 perf_log = perf_log, perf_state = perf_state,
                 analysis_cache = analysis_cache)

  spatialServer("spatial", datasets, active_id)
  exportServer("export", datasets, active_id)

  workflowServer("workflow", datasets, active_id,
                 preproc_specs = if (!is.null(preproc)) preproc$specs else reactive(list()))

  performanceServer("performance", datasets, active_id,
                    perf_log = perf_log, perf_state = perf_state,
                    analysis_cache = analysis_cache)

  helpServer("help", log_path = app_log_path)

  return(list(
    mse_data = mse_data
  ))
}
