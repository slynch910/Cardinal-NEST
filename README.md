# Cardinal NEST

A user-friendly graphical interface for mass spectrometry imaging (MSI) data analysis using the Cardinal 3 package.

## Overview

Cardinal NEST provides a visual workflow for MSI data analysis without requiring R programming expertise. It abstracts complex Cardinal functions into interactive controls, enabling end-to-end analysis:

- **Import**: Load imzML, Analyze 7.5, and other MSI data formats
- **Explore**: View metadata, spectra, and ion images
- **Preprocess**: Apply normalization, baseline correction, peak picking, and more
- **Analyze**: Perform PCA, classification, and spatial analysis
- **Export**: Save processed data and results

## Installation

```r
# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CardinalNEST")

# Or install development version
devtools::install_github("username/CardinalNEST")
```

## Quick Start

```r
library(CardinalNEST)

# Launch the application
runCardinalNEST()
```

## Features

### Data Management
- Import imzML and Analyze 7.5 files
- Interactive metadata explorer
- Dataset library with search and tagging
- Memory management for large datasets

### Visualization
- Interactive spectrum viewer
- Ion image viewer with multiple colormaps
- Metadata image mode: categorical pixelData maps (groups, clusters, segments)
- Figure Builder: collect and download spectra/ion images as a labelled grid
- ROI selection tools

### Preprocessing & Workflow Automation
- Drag-and-drop pipeline builder
- Save and load preprocessing workflows as JSON
- Export a workflow to a reproducible, headless R script
- Batch preprocessing: apply a saved workflow to multiple datasets in one run
- Normalization (TIC, RMS, reference)
- Baseline correction
- Peak picking and alignment
- Delayed processing for efficiency

### Statistical Analysis
- Univariate mean tests (`Cardinal::meansTest`) with FDR-aware feature ranking via `Cardinal::topFeatures()`
- Contrast testing (`Cardinal::contrastTest`) for post-hoc pairwise or custom contrasts on mean-test models, with configurable emmeans adjustment and feature-level FDR; results are sorted by the most significant contrast per feature
- PCA and NMF with scores, loadings/activation plots, variance/contribution, and spatial score maps
- Colocalization with a reference m/z or a selected region
- PLS-DA and OPLS-DA classification with scores scatter, spatial response image, and regression-coefficients plot, plus optional k-fold cross-validation (`Cardinal::crossValidate`)
- Supervised SSC (spatial shrunken centroids) classification against a class variable, with predicted-class map, class-probability image, and t-statistic feature ranking via `Cardinal::topFeatures()`, plus optional k-fold cross-validation; note that SSC feature ranking is by t-statistic (not FDR), unlike mean-test and contrast results
- k-means clustering
- Spatial segmentation (SSC and per-feature spatial-DGMM probabilistic mixtures) with feature ranking via `Cardinal::topFeatures()`
- Interactive ROI selection on pixel maps (rectangle or free-form polygon) with add/remove/clear-current-selection controls

### Performance & Parallel Computing
- Choose a BiocParallel backend (Serial, Snow, or Multicore where supported) and worker count, applied via `Cardinal::setCardinalBPPARAM()`
- Per-dataset memory panel showing pixel/feature counts, estimated in-memory size, and storage mode
- One-click "Convert to in-memory" for the active dataset to speed up repeated analysis
- Result cache for the Analyze tab: re-run identical inputs instantly
- Operation timing log for imports, preprocessing, and analysis runs

### Data Import
- imzML, Analyze 7.5, and CSV/tabular formats
- Multi-file upload: import several datasets in one action with auto-grouping of paired files
- Bulk/folder import: scan a server-side folder (recursively if desired) and import every detected dataset group with one click
- File-backed (shared memory) mode for large imzML datasets that exceed available RAM
- Configurable mass range and resolution filtering on import
- 10 GB maximum upload size for large datasets

## Documentation

- [User Guide](docs/user_guide/)
- [API Reference](docs/api_reference/)
- [Tutorials](docs/tutorials/)
- [Examples](docs/examples/)

## Development

Cardinal NEST is developed as a Cardinal 3 GUI application. All functionality leverages Cardinal's native capabilities rather than implementing general alternatives.

### Project Structure

```
cardinal-nest-shiny/
├── R/                      # R source code
│   ├── modules/           # Shiny modules
│   ├── utils/             # Utility functions
│   ├── ui/                # UI components
│   └── server/            # Server logic
├── docs/                   # Documentation
├── tests/                  # Test suite
└── inst/                   # Installed files
```

### Running Tests

The unit-test suite lives in `tests/unit/`. Launcher scripts resolve `Rscript`
automatically, so you do not need R on your `PATH`:

```powershell
# Windows (PowerShell) - full suite
powershell -ExecutionPolicy Bypass -File tests\run_tests.ps1

# Windows - a single file
powershell -ExecutionPolicy Bypass -File tests\run_tests.ps1 tests\unit\test_analysis.R
```

```bash
# Linux/macOS - full suite
bash tests/run_tests.sh
```

Resolution order for `Rscript`: the `RSCRIPT` environment variable, then `Rscript`
on the `PATH`, then the newest `R-*` install under `Program Files` (Windows) or the
usual system locations (Linux/macOS). To force a specific interpreter, set
`RSCRIPT` to its full path before running. The launcher exits `0` on success or
with the number of failures otherwise.

## Contributing

Contributions are welcome! Please read our contributing guidelines and code standards in `docs/CODE_STANDARDS.md`.

## License

Artistic-2.0

## Citation

If you use Cardinal NEST in your research, please cite:

> Cardinal NEST: A Graphical User Interface for Cardinal 3. [Year].

## Acknowledgments

Cardinal NEST is built on top of the Cardinal package and the broader Bioconductor ecosystem.
