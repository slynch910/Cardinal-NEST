# Changelog

All notable changes to Cardinal NEST will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- README rewritten for the `v0.2.0-preview` public release: preview disclaimer, `renv::activate()` + `renv::restore()` install from `renv.lock`, local launch, remote/cluster launch via SSH tunnel (direct host and compute-node-behind-login-node), a Known Issues list, and GitHub Issues as the feedback channel; removed the non-existent Bioconductor/`runCardinalNEST()` install instructions and dead `docs/` links

## [0.2.0] - 2026-10-14

### Changed
- Normalized source comments: removed full-width `#---`/`#===` section banners, removed restating boilerplate (`# Module UI/server function`), dropped "Mirrors X" cross-references, and condensed verbose header/comment blocks across the `R/` modules and `app.R`
- Spatial tab help text now uses the user-facing label "Spatial Shrunken Centroids (SSC)" instead of the package function name `Cardinal::spatialShrunkenCentroids`
- Feature ranking in `run_means_test`, `ssc_top_features`, and `segmentation_top_features` now uses `Cardinal::topFeatures()` for consistent ranking; note that SSC/segmentation tables are ranked by the shrunken t-statistic (no FDR), while mean-test and contrast tables are FDR-aware
- `run_contrast_test` now sorts features by the minimum p-value across all contrasts, not just the first contrast
- Cross-validation (`run_cross_validation`) now uses reproducible, shuffled, stratified fold assignments and canonicalizes PLS/SSC score matrices to a metrics-as-rows layout
- Cross-validation score metric renamed from `accuracy` to `mean_recall` to reflect that the reported value is the per-class recall averaged across folds
- `run_cross_validation` defaults for `r` and `s` now match `run_ssc_classification` (`r = 2`, `s = 9`)
- `read_imzml` now validates `mass_range` and `resolution` before passing them to `Cardinal::readMSIData`
- `build_imzml_opts()` helper centralizes Shiny input normalization for imzML import settings

### Fixed
- Export tab no longer drops added `pixelData` metadata columns (cluster/segment/ROI/PCA-NMF-DA score columns such as `kmeans_cluster`, `spatial_segment`, `dgmm_segment`, `ssc_class`, `roi`, `PC1...`, `NMF_...`): imzML, Analyze 7.5, and CSV exports now write a per-pixel `*.metadata.csv` sidecar (bundled into the download zip) so the metadata persists in the exported file even without re-import; rds already carried it inside the serialized object. New pure helpers `added_metadata_columns()`, `pixel_metadata_frame()`, `write_metadata_sidecar()`, and `export_delivers_bundle()`, with unit tests in `tests/unit/test_export.R`
- Coerce `r` and `s` to `numeric` before passing them to `Cardinal::spatialShrunkenCentroids`, preventing the `REAL() can only be applied to a 'numeric', not a 'integer'` runtime error
- Fixed analysis engines (PCA, NMF, PLS-DA, OPLS-DA, SSC, k-means, cross-validation) on file-backed datasets by coercing spectra to an in-memory matrix before `run_analysis`
- Fixed cross-validation `cv@model` extraction so the result list remains safely usable with `[[` despite the `SpatialCV` object's custom `$` method

### Added
- `renv.lock` pinning R, Cardinal, Bioconductor, and all other package dependencies so installs can be reproduced with `renv::restore()`
- Colon range syntax for index inputs in the Visualize tab: Spectrum Index and Feature Index now accept inclusive `a:b` ranges mixed with commas (e.g. `1:10, 30`). Reversed, incomplete, or non-whole ranges are rejected with a clear message and duplicates are removed while keeping order. m/z input stays comma-only and rejects ranges. `parse_numeric_list()` moved to `R/utils/helpers.R` as a pure helper, with unit tests in `tests/unit/test_visualization.R`
- R data (.rds) import in the Import tab: load a Cardinal MSI object previously saved via the Export tab (or `saveRDS`), completing the RDS export/import round-trip; supports single/multiple upload and folder bulk import, with `read_rds_msi` pure helper that validates the object is an MSI type and coerces `MSImagingArrays` to `MSImagingExperiment`
- Workflow tab (`R/modules/workflow/workflow_module.R`): save/load preprocessing pipelines as versioned JSON, export a reproducible R script for headless execution, and batch-apply a workflow across multiple datasets with a per-dataset status table
- Bulk/folder import in the Import tab: scan a server-side folder for imzML, Analyze 7.5, or CSV files and import every detected dataset group with one click, plus an optional recursive subfolder scan
- `jsonlite` dependency for workflow serialization
- Performance tab (`R/modules/performance/performance_module.R`): configure a BiocParallel backend (Serial/Snow/Multicore) with worker count, view estimated memory/storage mode for each dataset, convert the active dataset to in-memory with one click, enable/clear an analysis result cache for instant re-runs, and monitor import/preprocess/analysis timings in a shared log

### Added
- Contrast testing in the Analyze tab: post-hoc pairwise (or custom) contrasts on `Cardinal::meansTest` models via `Cardinal::contrastTest`, with configurable contrast method, emmeans p-value adjustment, and feature-level FDR correction, plus a significance plot
- Cross-validation for PLS-DA and SSC in the Analyze tab: optional k-fold cross-validation via `Cardinal::crossValidate`, with per-class recall/precision bar chart and mean recall reported in the status message
- `read_imzml` pure helper for testable file-backed (shared) or in-memory imzML import with optional mass range and resolution filtering
- Raised Shiny max upload size to 10 GB (`options(shiny.maxRequestSize)`) to accommodate large imzML/ibd datasets
- Supervised SSC classification in the Analyze tab: `Cardinal::spatialShrunkenCentroids` trained with `y = class` from the existing grouping infrastructure, with a "Plot view" selector for predicted-class map, class-probability image (`image(model, type='probability')`), and shrunken t-statistic plot (`plot(model, type='statistic')`), plus a save-to-dataset flow that stores the predicted class as `ssc_class`
- PLS-DA and OPLS-DA outputs in the Analyze tab: the existing discriminant-analysis results now support a "Plot view" selector to switch between scores scatter, spatial response image (`image(model, type='response')`), and regression-coefficients plot (`plot(model, type='coefficients')`)
- Spatial-DGMM probabilistic segmentation in the Spatial tab: run `Cardinal::spatialDGMM` on a chosen feature (m/z), visualize its per-pixel Gaussian-mixture segment map, and save the segment assignment as a new `pixelData()` column (`dgmm_segment`)
- "Clear current selection" button in the Spatial tab ROI selection: removes the in-progress rectangle brush and polygon outline from the map without discarding the committed ROI; switching shape modes now also resets the brush and outline
- Metadata image mode in the Visualize tab: render any categorical `pixelData()` column (group labels, clusters, segments, ROI membership) as a spatial map with each pixel colored by its category
- Free-form (polygon) ROI selection in the Spatial tab: click vertices on the pixel map to trace an irregular outline, then add or remove the enclosed pixels from the ROI, with undo/clear outline controls
- Colocalization tool in the Analyze tab: rank features co-localized with a reference m/z (via `Cardinal::colocalized`) or with a selected class-level region (via feature-image correlation against the region mask), and display the top co-localized ions as ion images
- NMF (non-negative matrix factorization) analysis in the Analyze tab via `Cardinal::NMF`, with scores scatter plots, activation (loadings) vs m/z, per-component contribution bar chart, and spatial score maps (`image(nmf, type = 'x')`); NMF scores can be saved as new `pixelData` columns on a new in-memory dataset
- PCA completion: added PCA loadings plots (`plot_pca_loadings()`), variance-explained bar chart (`plot_pca_variance()`), and spatial score maps (`image(pca, type = 'x')`); PCA and NMF plots are now switchable via a "Plot view" selector (scores / loadings / component contribution / spatial score map)
- ROI selection in the Spatial tab: interactively drag rectangles on the pixel map to add/remove pixels from a named region of interest, view the live selection count, export ROI membership to CSV, and save the ROI as a new pixelData column on a new in-memory dataset
- `pixels_in_brush()`, `make_roi_labels()`, `roi_pixel_table()`, `plot_roi()`, and `make_roi_record()` pure helpers for ROI selection, with unit tests covering brushing, membership labelling, plotting, record creation, and the add/remove/clear/save server wiring (`tests/unit/test_spatial.R`)
- Figure Builder in the Visualize tab: collect displayed spectra/ion-image plots into a grid, add/remove/clear panels, preview, and download the assembled figure as PNG or PDF with optional per-panel labels "(a)", "(b)", ...
- Figure Builder panels are individually captured as grobs (`gridGraphics::echoGrob`) and composed with `gridExtra::arrangeGrob`, which lets spectra and ion images with different coordinate systems share one figure without breaking
- Added per-panel editable labels and an explicit "Preview" modal so the full assembled figure is shown only on demand, not as a real-time inline preview
- Added `split_display_plot()` to itemize multi-panel `vizi_facets` (e.g., several non-superposed spectra) into separate, individually labelled collection entries
- Unit tests for the Figure Builder helpers and add/remove/clear/preview/download server wiring (`tests/unit/test_visualization.R`)
- Spatial Analysis tab (`R/modules/spatial/spatial_module.R`): a dedicated tab for spatial/coordinate-based analysis, separate from the statistical `Analyze` tab
- Spatial segmentation via `Cardinal::spatialShrunkenCentroids` with user-configurable `k`, `r`, `s`, and Gaussian/adaptive spatial weights; displays a segment map, a ranked top-features table, and supports saving per-pixel segments as a new in-memory dataset
- `run_spatial_segmentation()`, `segmentation_top_features()`, `segmentation_pixel_table()`, `plot_segmentation()`, and `make_segmentation_record()` pure helpers for the Spatial module
- Unit tests for the Spatial module helpers and server wiring (`tests/unit/test_spatial.R`)
- Help tab (`R/modules/help/help_module.R`) with a `Log` panel that captures R console `output` and `message`/`warning` text via `sink()` and streams it back into the app with `reactivePoll()`; `Vignette` panel is a placeholder for future user guides
- Supervised discriminant analysis in the Analyze tab: PLS-DA (`Cardinal::PLS`) and OPLS-DA (`Cardinal::OPLS`) built on the existing grouping infrastructure (`run_discriminant_analysis()`)
- Discriminant results: latent-variable scores/loadings, per-pixel predicted class, training accuracy, and confusion matrix; interactive component pickers for the scores plot (`plot_da_scores()`)
- Save PLS-DA/OPLS-DA scores and predicted class back onto the dataset as a new in-memory dataset (`make_da_record()`)
- Unsupervised clustering in the Analyze tab: spatially-aware k-means (`Cardinal::spatialKMeans`)
- Cluster results: per-pixel cluster labels, cluster sizes, and a spatial cluster map (`plot_clustering()`); cluster labels can be saved back onto the dataset as a new in-memory dataset (`make_clustering_record()`)
- Unit tests for PLS-DA/OPLS-DA helpers, plotting, score-record creation, `testServer()` run/save wiring, and for k-means clustering helpers and save flow (`tests/unit/test_analysis.R`)
- Full CSV/tabular import: parse `x`/`y`/`run` columns and numeric m/z headers into an `MSImagingExperiment` (`parse_csv_msi()`, `build_msi_from_csv()`, `read_csv_msi()`)
- Multi-dataset import: uploading multiple imzML/Analyze pairs or several tabular files imports each as a separate dataset in one action (`import_group_files()`)
- Progress indicators (`withProgress`/`incProgress`) for imports and for univariate/PCA analysis runs
- Import history log tracking every import attempt (time, name, format, status, dimensions, details) in a session table
- Data preview card showing the active dataset's mean spectrum and first pixel coordinates
- `summarize_msi()` helper to build the standard dataset summary
- Unit tests for CSV parsing, multi-file grouping, summary generation, and CSV round-trip (`tests/unit/test_data_import.R`)
- Export module for saving the active dataset to disk in `rds`, `imzML`, `img` (Analyze 7.5), and `csv` formats (`R/modules/export/export_module.R`)
- `export_basename()`, `export_csv_matrix()`, `export_dataset()`, and `write_export_download()` pure helpers for filename sanitization and multi-format disk export
- Integration of `exportUI()`/`exportServer()` into `app.R`, `R/ui/main_ui.R`, and `R/server/main_server.R`
- Export unit tests covering format metadata, filename sanitization, CSV matrix shape, `imzML`/`img`/`rds` file output, `.zip` packaging, `rds` round-trip, and `testServer()` preview/active-dataset wiring (`tests/unit/test_export.R`)
- Visualization module for plotting spectra and rendering ion images (`R/modules/visualization/visualization_module.R`)
- Spectra plotting by index, pixel coordinate, or summary (`Cardinal::plot`)
- Ion image rendering by m/z, feature index, or pixel summary (`Cardinal::image`, `Cardinal::summarizePixels`)
- Multi-value input support for spectrum indices and ion m/z / feature indices (comma-separated)
- Superpose controls for overlaying multiple spectra or ion images
- Summary spectra using `Cardinal::summarizeFeatures` and pixel summary ion images (TIC/Mean)
- `numericRangeInput` m/z range selector in spectra viewer
- Plot loading spinner (`shinycssloaders`) and responsive plot sizing
- Comprehensive unit tests for the visualization module (`tests/unit/test_visualization.R`)
- Basic Shiny app structure with main_ui.R and main_server.R
- Tabbed navigation interface (Import, Explore, Preprocess, Visualize, Analyze, Export)
- Reactive data storage for MSImagingExperiment objects
- Helper functions placeholder file
- Data import module with file selection, format options, and advanced settings
- File validation and error handling for imzML import
- Data summary display showing pixels, features, m/z range, and memory mode
- Integration with Cardinal's readMSIData function for file import
- Unit tests for data import module using testthat
- Test coverage for file validation, parameter handling, and data structure verification
- File pair validation tests for imzML and analyze formats
- Memory parameter handling tests
- Initial project structure and framework
- Project planning documentation with Cardinal 3 principles
- Development standards for unit testing and change logging
- Infrastructure consistency guidelines for S3/S4 class usage

### Changed
- Upgraded the app to Bootstrap 5 using `bslib::bs_theme(version = 5)` in `navbarPage`
- Replaced `wellPanel()` with `bslib::card()` across data import, metadata, and visualization modules
- Refactored visualization module with helper functions for input parsing, shared plot rendering, and axis expansion
- Renamed visualization inputs for consistency: `plot_mode`, `spec_*`, `viz_result`, `spec_show_key`
- Moved package loads for `Cardinal`, `bslib`, `shinyjs`, `shinyWidgets`, and `shinycssloaders` to `app.R`
- Improved UI responsiveness by switching from sidebarLayout to navbarPage for better window resizing
- Changed from single tab_content output to individual tab-specific outputs for better structure
- Integrated data import module into main application
- Fixed file upload to accept all file types (removed restrictive accept parameter)
- Added dynamic file input label that changes based on selected format
- Fixed deprecated parameter by replacing attach.only with memory parameter for readMSIData
- Added format-specific file validation for Analyze 7.5 and CSV formats
- Restricted advanced settings to only show for imzML format
- Updated app.R to clarify it's a development script for package structure
- Implemented multi-file upload support for paired files (imzML+ibd, hdr+img)
- Added file pair validation logic to ensure both required files are uploaded together
- Implemented temporary directory management with proper file naming for Cardinal compatibility
- Added automatic cleanup of temporary directories after import or on error
- Added centroided status to data summary display for user information (not tested as Cardinal doesn't require centroided data)
- Added comprehensive code breakdown documentation with function explanations and documentation links
- Removed redundant data_import_server.R file (functionality consolidated into data_import_module.R)
- Implemented metadata management module with pixel and feature data exploration
- Added interactive data tables for pixelData() and featureData() with filtering capabilities
- Added spatial coordinate filtering (x, y ranges) for pixel metadata
- Added m/z range filtering for feature metadata
- Added memory management controls for in-memory vs file-backed data
- Integrated metadata module into Explore tab
- Created unit tests for metadata management module
- Fixed metadata filter logic to support min-only, max-only, or both bounds for x, y, and m/z ranges
- Fixed metadata filters to only apply when the Apply button is pressed (no longer triggers on enter key or in real-time)
- Implemented multi-dataset management: imported datasets are now stored together instead of only keeping the most recent
- Added a sortable/filterable dataset table in the Import tab to browse, select, and delete imported datasets
- Added active dataset selection so any imported dataset can be made active for exploration
- Added active dataset indicator to the Explore tab
- Refactored data import and metadata modules to share a central dataset store (datasets + active_id)
- Created unit tests for dataset management (store add/delete, id generation, overview building, active reassignment)
- Fixed flawed metadata unit tests: coord() returns an S4 DataFrame (not base matrix/data.frame), and m/z filter test now uses independent min/max bounds matching the module logic
- Updated IMPLEMENTATION.md and CODE_BREAKDOWN.md to document metadata and dataset management systems

### Fixed
- Fixed spectrum y-axis and x-axis ticks so they extend up to the data maximum by expanding plot limits to `pretty()` bounds
- Fixed `numericInput` step defaults so tolerance, m/z range, and line width accept decimal values
- Wrapped `Cardinal::plot` and `Cardinal::image` calls in `tryCatch` so invalid user inputs surface messages instead of crashing the app
- Changed `btn-default` to `btn-secondary` in metadata module for Bootstrap 5 compatibility
- Fixed Analyze 7.5 import failure by validating that .hdr, .img, and .t2m files are all uploaded together; Cardinal's readMSIData() requires all three files for Analyze datasets
- Updated unit tests and documentation to reflect the Analyze 7.5 .t2m requirement

## [0.1.0] - 2026-07-14

### Added
- Project initialization
- Basic directory structure
- Package DESCRIPTION and NAMESPACE files
- README.md with project overview
