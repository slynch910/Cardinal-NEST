# Cardinal NEST - main server

main_server <- function(input, output, session,
                        log_path = file.path(getwd(), "logs", "app.log"),
                        vignette_dir = tutorials_dir()) {

  # Dataset store: named list of records (metadata + MSImagingExperiment).
  datasets <- reactiveVal(list())
  active_id <- reactiveVal(NULL)

  perf_log <- reactiveVal(NULL)
  perf_state <- reactiveValues(cache_enabled = TRUE, global_seed = default_seed())
  # Reaches modules with a fixed signature (spatialServer) via userData.
  session$userData$perf_state <- perf_state
  analysis_cache <- new_result_cache()

  # Session operation log (Help tab); published for module logging helpers.
  op_log <- reactiveVal(list())
  session$userData$op_log <- op_log

  # Log create/remove/rename diffs of the dataset store.
  store_snapshot_rv <- reactiveVal(list())
  observe({
    ds <- datasets()
    for (ev in dataset_store_events(store_snapshot_rv(), ds)) {
      log_session_event(ev)
    }
    store_snapshot_rv(store_snapshot(ds))
  })

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

  helpServer("help", log_path = log_path, vignette_dir = vignette_dir,
             op_log = op_log)

  return(list(
    mse_data = mse_data
  ))
}
