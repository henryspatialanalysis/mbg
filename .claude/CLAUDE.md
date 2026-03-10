# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install package locally
make install        # runs devtools::install()

# Rebuild roxygen2 documentation and NAMESPACE
make build-docs     # runs devtools::document()

# R CMD CHECK equivalent (standard R package check)
Rscript -e "devtools::check()"

# Run a single function's tests interactively
Rscript -e "devtools::load_all(); <function_call>"

# Rebuild pkgdown site
make deploy-docs    # builds and deploys to GitHub Pages
```

There are no formal unit tests (no `tests/` directory). Validation is done through vignettes in `vignettes/`.

## Architecture

**mbg** is an R package for Model-Based Geostatistics — estimating continuous spatial surfaces from point-referenced observations combined with raster covariates.

### Core Components

**`MbgModelRunner` R6 class** ([R/MbgModelRunner.R](R/MbgModelRunner.R)) is the main orchestrator. It holds all inputs as fields and calls the modular functions below in sequence. All functions can also be called standalone.

**Pipeline stages:**

1. **Input loading** — `load_covariates()`, `build_id_raster()`, `build_aggregation_table()` prepare rasters and pixel-to-polygon mappings.

2. **ML stacking (optional)** — `run_regression_submodels()` ([R/run_ml_models.R](R/run_ml_models.R)) fits caret models (Random Forest, GAM, etc.) to produce predictive raster surfaces used as additional covariates in the INLA model.

3. **INLA preparation** — `prepare_inla_data_stack()` ([R/prepare_inla_data_stack.R](R/prepare_inla_data_stack.R)) builds the INLA data stack with fixed effects (covariates), SPDE spatial mesh (Gaussian process), optional admin-level IID random effects, and optional nugget term.

4. **Model fitting** — `fit_inla_model()` ([R/fit_inla_model.R](R/fit_inla_model.R)) wraps `INLA::inla()` with PC priors and sensible defaults.

5. **Prediction** — `generate_cell_draws_and_summarize()` ([R/generate_cell_draws_and_summarize.R](R/generate_cell_draws_and_summarize.R)) draws 250 posterior predictive samples per pixel, applies inverse link function, and produces mean/median/UI rasters.

6. **Aggregation** — `aggregate_draws_to_polygons()` ([R/aggregate_grid_to_polygons.R](R/aggregate_grid_to_polygons.R)) propagates uncertainty through aggregation from raster cells to polygon-level summaries.

### Key Design Decisions

- **INLA is a Suggested dependency** (not required), checked at runtime in [R/zzz.R](R/zzz.R). The package can be installed without INLA.
- **terra >= 1.9.1 required** for efficient aggregation table support.
- Raster operations use **terra**, vector/polygon operations use **sf**, tabular data uses **data.table**.
- The `id_raster` (created by `build_id_raster()`) defines the spatial template — all rasters and predictions align to it.
- Aggregation uses fractional pixel allocation, preserving uncertainty through the `cell_draws` matrix (pixels × samples).

### Example Data

`inst/extdata/` contains Benin communes boundaries and child stunting observations, used in the vignettes as a worked example.
