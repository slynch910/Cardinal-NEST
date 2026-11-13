# Standalone parallel backend check (no Shiny session needed).
# Confirms a BiocParallel backend uses separate worker processes and times
# Cardinal::summarizeFeatures() serially versus on that backend, using the
# same helpers as the Performance tab's "Test backend" button.
#
# Usage, from the repository folder (the one containing app.R):
#   Rscript inst/scripts/verify_backend.R backend=snow workers=4
#   Rscript inst/scripts/verify_backend.R backend=multicore workers=8 data=/path/run.imzML
# Arguments (all optional):
#   backend  serial | snow | multicore (default snow; multicore is not available on Windows)
#   workers  number of workers (default: parallel::detectCores() - 1, at least 2)
#   seed     RNG seed for the backend (default: the global seed, normally 42)
#   data     imzML file to time; without it a small simulated image is used
# Exits with status 1 when a parallel backend did not use separate workers.

suppressPackageStartupMessages({
  library(Cardinal)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
opts <- list()
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2L) opts[[kv[1]]] <- kv[2]
}

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."))
} else {
  getwd()
}
source(file.path(root, "R", "helpers.R"))
source(file.path(root, "R", "performance_module.R"))

backend <- opts$backend %||% "snow"
workers <- opts$workers %||% max(2L, parallel::detectCores() - 1L)
seed <- resolve_seed(opts$seed)
param <- make_bpparam(backend, workers, seed = seed)

data <- if (!is.null(opts$data)) {
  readMSIData(opts$data, verbose = FALSE)
} else {
  simulateImage(preset = 1L, dim = c(64L, 64L), verbose = FALSE,
                BPPARAM = SerialParam(RNGseed = seed))
}

cat(sprintf("R %s, Cardinal %s, BiocParallel %s\n", getRversion(),
            packageVersion("Cardinal"), packageVersion("BiocParallel")))
cat(sprintf("Data: %d features x %d pixels (%s)\n", nrow(data), ncol(data),
            if (is.null(opts$data)) "simulated" else opts$data))
cat(sprintf("Main process PID: %d\n", Sys.getpid()))

res <- test_backend(param, data)
summary <- backend_test_summary(res)
cat(sprintf("[%s] %s\n", toupper(summary$type), summary$message))

if (identical(summary$type, "warning")) quit(status = 1L)
