# Package index

## MBG Model Runner

Model runner class

- [`MbgModelRunner`](https://henryspatialanalysis.github.io/mbg/reference/MbgModelRunner.md)
  : MBG model runner class

## MBG core inputs

Functions to build the ID raster and load model covariates

- [`build_id_raster()`](https://henryspatialanalysis.github.io/mbg/reference/build_id_raster.md)
  : Build ID raster
- [`load_covariates()`](https://henryspatialanalysis.github.io/mbg/reference/load_covariates.md)
  : Load covariates
- [`make_world_template_raster()`](https://henryspatialanalysis.github.io/mbg/reference/make_world_template_raster.md)
  : Make world template raster

## MBG model fitting and prediction

Individual functions for model fitting and prediction. These can also be
executed through the `MbgModelRunner` class.

- [`fit_inla_model()`](https://henryspatialanalysis.github.io/mbg/reference/fit_inla_model.md)
  : Fit INLA model
- [`prepare_inla_data_stack()`](https://henryspatialanalysis.github.io/mbg/reference/prepare_inla_data_stack.md)
  : Prepare data stack for INLA
- [`run_regression_submodels()`](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md)
  : Run regression sub-models
- [`generate_cell_draws_and_summarize()`](https://henryspatialanalysis.github.io/mbg/reference/generate_cell_draws_and_summarize.md)
  : Generate cell draws and summary rasters from INLA model
- [`summarize_draws()`](https://henryspatialanalysis.github.io/mbg/reference/summarize_draws.md)
  : Summarize draws
- [`values_to_raster()`](https://henryspatialanalysis.github.io/mbg/reference/values_to_raster.md)
  : Insert values into a raster

## Model validation

Functions for generating predictive validity metrics. These can also be
executed through the `MbgModelRunner` class.

- [`log_posterior_density()`](https://henryspatialanalysis.github.io/mbg/reference/log_posterior_density.md)
  : Generate log posterior predictive density from a geostatistical
  surface onto point data
- [`rmse_raster_to_point()`](https://henryspatialanalysis.github.io/mbg/reference/rmse_raster_to_point.md)
  : Generate RMSE from an estimated raster surface and point data

## Pixel-to-polygon aggregation

Functions to aggregate rasters and pixel-level predictive draw matrices
to polygon boundaries. These can also be executed through the
`MbgModelRunner` class.

- [`aggregate_draws_to_polygons()`](https://henryspatialanalysis.github.io/mbg/reference/aggregate_draws_to_polygons.md)
  : Aggregate grid cell draws to polygons
- [`aggregate_raster_to_polygons()`](https://henryspatialanalysis.github.io/mbg/reference/aggregate_raster_to_polygons.md)
  : Aggregate a raster to polygons
- [`build_aggregation_table()`](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md)
  : Build aggregation table

## Other

Miscellaneous helper functions

- [`dissolve_sf_by_attribute()`](https://henryspatialanalysis.github.io/mbg/reference/dissolve_sf_by_attribute.md)
  : Dissolve sf object by attribute
- [`.onAttach()`](https://henryspatialanalysis.github.io/mbg/reference/dot-onAttach.md)
  : Behavior when attaching the mbg package
- [`logging_get_timer_log()`](https://henryspatialanalysis.github.io/mbg/reference/logging_get_timer_log.md)
  : Get timer log
- [`logging_start_timer()`](https://henryspatialanalysis.github.io/mbg/reference/logging_start_timer.md)
  : Start logging timer
- [`logging_stop_timer()`](https://henryspatialanalysis.github.io/mbg/reference/logging_stop_timer.md)
  : End logging timer
- [`make_time_stamp()`](https://henryspatialanalysis.github.io/mbg/reference/make_time_stamp.md)
  : Make time stamp
- [`vif_covariate_select()`](https://henryspatialanalysis.github.io/mbg/reference/vif_covariate_select.md)
  : Run variance inflation factor (VIF) selection on input covariates
