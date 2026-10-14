# Cardinal NEST

A user-friendly graphical interface for mass spectrometry imaging (MSI) data analysis using the Cardinal 3 package.

> **Preview release (`v0.2.0-preview`)**: This is an early build meant to let lab members try the app. Expect rough edges, and expect behavior and interfaces to change before the stable release. Please check [Known Issues](#known-issues) and report problems through [GitHub Issues](https://github.com/slynch910/Cardinal-NEST/issues).

## Overview

Cardinal NEST provides a visual workflow for MSI data analysis without requiring R programming expertise. It abstracts complex Cardinal functions into interactive controls, enabling end-to-end analysis:

- **Import**: Load imzML, Analyze 7.5, and other MSI data formats
- **Explore**: View metadata, spectra, and ion images
- **Preprocess**: Apply normalization, baseline correction, peak picking, and more
- **Analyze**: Perform PCA, classification, and spatial analysis
- **Export**: Save processed data and results

## Installation

The preview runs from a cloned copy of the repository. Package versions are pinned in `renv.lock` (R 4.6, Bioconductor 3.23, Cardinal 3.14), so `renv` gives every machine the same environment.

### Requirements

- **R 4.6.x**. `renv` warns if your R version differs from the lockfile.
- **Git**, or download the repository as a ZIP.
- **Build tools** for packages compiled from source: Rtools on Windows, Xcode Command Line Tools on macOS, or a C/C++/Fortran toolchain on Linux and clusters (for example `module load gcc`).

### Steps

```bash
git clone https://github.com/slynch910/Cardinal-NEST.git
cd Cardinal-NEST
```

Then, from an R session started **in the `Cardinal-NEST` folder**:

```r
install.packages("renv")
renv::activate()  # creates a project-private library for this folder
```

**Restart R** (or start a new R session in the `Cardinal-NEST` folder), then:

```r
renv::restore(prompt = FALSE)  # installs the exact versions in renv.lock
```

The restart matters: if `restore()` runs in the same session as `activate()`, the project library is not yet active and restore can incorrectly report "already synchronized" while installing nothing.

The first restore compiles some packages and can take a while. Packages go into the project library, not your personal R library. Any R session started in this folder uses that library automatically.

## Launching

Always launch from the repository folder (the one that contains `app.R`).

### Local machine

```bash
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
```

Or, inside R: `shiny::runApp(".")`. By default the app listens on `http://127.0.0.1:8080`.

### Remote server or cluster (SSH tunnel)

Run the app on the remote machine and view it in your local browser through an SSH tunnel. Pick a port nobody else is using (the examples use `8765`).

**Case 1: the app runs on the machine you SSH into** (a workstation, homelab, or login node where running jobs is allowed).

On the remote machine:

```bash
cd Cardinal-NEST
Rscript -e "shiny::runApp('.', host = '127.0.0.1', port = 8765, launch.browser = FALSE)"
```

On your local machine:

```bash
ssh -N -L 8765:127.0.0.1:8765 user@remote-host
```

Then open `http://localhost:8765` in your local browser.

**Case 2: the app runs on a compute node behind a login node** (for example an interactive Slurm job).

On the compute node (inside your interactive job):

```bash
hostname              # note the node name, e.g. node042
cd Cardinal-NEST
Rscript -e "shiny::runApp('.', host = '0.0.0.0', port = 8765, launch.browser = FALSE)"
```

On your local machine, tunnel through the login node to the compute node:

```bash
ssh -N -L 8765:node042:8765 user@login-node
```

Then open `http://localhost:8765`.

Notes:

- `host = '0.0.0.0'` makes the app reachable from other machines on the cluster network, and the app has **no authentication**. Use it only on trusted internal networks, and prefer `127.0.0.1` (Case 1) whenever possible.
- Keep the remote session alive while you work, for example with `tmux` or `screen`, or inside the batch job.
- Logs are written to `logs/app.log` in the repository folder and can be viewed in the **Help** tab.
- On a cluster, files that are already there can be loaded with the server-side folder import on the **Import** tab. This avoids uploading them through the browser.

## Known Issues

These are known limitations of the preview. They are planned for upcoming releases:

- **File-backed datasets become in-memory** after analysis, spatial, and export operations, so very large datasets can use a lot of RAM.
- **Randomized analyses are not seeded yet** (PCA, NMF, spatial k-means, spatial-DGMM), so repeated runs can differ slightly.
- **Batch preprocessing results are kept only for the session** and are not written to disk automatically. Use the **Export** tab to save them.
- **The Figure Builder uses a fixed panel width**, so panels with mixed aspect ratios can look stretched.
- **No package install yet**: `runCardinalNEST()` and `remotes::install_github()` installation are planned for the stable release. For now, launch from the cloned folder as shown above.
- **No example data ships with the repository**. Use your own imzML, Analyze 7.5, or CSV data.

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

- [Tutorials](docs/tutorials/)

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

The test suite is available only in the development repository. The public release does not include `tests/`. The unit-test suite lives in `tests/unit/`. Launcher scripts resolve `Rscript`
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

Bug reports and feature requests are welcome through [GitHub Issues](https://github.com/slynch910/Cardinal-NEST/issues). Please include your operating system, R version, and the relevant part of `logs/app.log`.

## License

Artistic-2.0

## Citation

If you use Cardinal NEST in your research, please cite:

> Cardinal NEST: A Graphical User Interface for Cardinal 3. [Year].

## Acknowledgments

Cardinal NEST is built on top of the Cardinal package and the broader Bioconductor ecosystem.
