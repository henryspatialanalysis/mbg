# MBG model runner class

R6 class to run a full MBG model and make predictions.

## Details

To see examples of this object, run
[`vignette('mbg')`](https://henryspatialanalysis.github.io/mbg/articles/mbg.md)

## Public fields

- `input_data`:

  ([data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html))  
  Table containing at least the following fields:  

  - x (`numeric`) location longitude in decimal degrees  

  - y (`numeric`) location latitude in decimal degrees  

  - indicator (`integer`) The number of events in the population  

  - samplesize (`integer`) The total population, denominator for
    `indicator`

- `id_raster`:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))  
  raster showing the total area that will be predicted using this model.

- `covariate_rasters`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A list containing all predictor covariates. Each covariate is a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the same extent and dimensions as `id_raster`.

- `aggregation_table`:

  ([data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html))  
  A table created by
  [build_aggregation_table](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md),
  used to link each grid cell to higher-level administrative units.

- `aggregation_levels`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list: for each named item, the name is the label for that
  aggregation level, and the value is a character vector of all fields
  in the original polygons to be used for aggregation at that level.

- `population_raster`:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))  
  A raster giving population for each grid cell, to be used for
  population-weighted aggregation from grid cells to polygon boundaries.
  Should have the same dimensions as `id_raster`. If no population
  raster is passed and the results are aggregated, aggregation will be
  by simple mean rather than population-weighted mean

- `admin_bounds`:

  ([sf::sf](https://r-spatial.github.io/sf/reference/sf.html))  
  Polygons showing the boundaries of administrative divisions within the
  study region. Only required if `use_admin_effect` OR
  `stacking_use_admin_bounds` is `TRUE`.

- `admin_bounds_id`:

  (`character`)  
  Field containing unique identifiers for `admin_bounds`, if passed.

- `use_covariates`:

  (`logical(1)`)  
  Should covariate effects be included in the predictive model?

- `use_gp`:

  (`logical(1)`)  
  Should a smoothed spatial surface be included in the predictive model?

- `use_admin_effect`:

  (`logical(1)`)  
  Should IID administrative-level effects be included in the predictive
  model?

- `use_nugget`:

  (`logical(1)`)  
  Should an IID effect by pixel be included in the predictive model?

- `use_stacking`:

  (`logical(1)`)  
  Should machine learning submodels be trained to relate the covariate
  rasters with the outcome data? Only run if `use_covariates` is `TRUE`.

- `stacking_model_settings`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list of submodels to be run. For more information about this
  term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `stacking_cv_settings`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  How should the stacking submodels be cross-validated? For more
  information about this term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is `TRUE`.

- `stacking_use_admin_bounds`:

  (`logical(1)`)  
  Should admin boundaries be included as features in the stacking
  submodels? For more information about this term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `stacking_prediction_range`:

  (`logical(1)`)  
  Range of possible predictions for the stacking submodels. For more
  information about this term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `mesh_max_edge`:

  (`numeric(2)` or NULL)  
  Maximum size of the INLA SPDE mesh inside (1) and outside (2) of the
  modeled region. Only considered if `use_gp` is TRUE.

- `mesh_cutoff`:

  (`numeric(1)`)  
  Minimum size of the INLA mesh, usually reached in data-dense areas.
  Only considered if `use_gp` is TRUE.

- `spde_integrate_to_zero`:

  (`boolean(1)`)  
  Should the 'volume' under the SPDE mesh integrate to zero? Only
  considered if `use_gp` is TRUE.

- `prior_spde_range`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list specifying the penalized complexity prior for the SPDE
  range. The two named items are "threshold", the test threshold (set as
  a proportion of the overall mesh extent), and "prob_below", the prior
  probability that the value is BELOW that range threshold. The function
  automatically converts "threshold" from a proportion of the overall
  mesh extent into a distance. Only considered if `use_gp` is TRUE.

- `prior_spde_sigma`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list specifying the penalized complexity prior for sigma
  (standard deviation) of the SPDE object. The two named items are
  "threshold", the test threshold for the standard deviation, and
  "prob_above", the prior probability that sigma will EXCEED that
  threshold. Only considered if `use_gp` is TRUE

- `prior_nugget`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list specifying the penalized complexity prior for the nugget
  term. The two named items are "threshold", the test threshold for the
  nugget standard deviation, and "prob_above", the prior probability
  that the standard deviation will EXCEED that threshold. Only
  considered if `use_nugget` is TRUE.

- `prior_admin_effect`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list specifying the penalized complexity prior for the
  admin-level IID term. The two named items are "threshold", the test
  threshold for the standard deviation of admin-level effects, and
  "prob_above", the prior probability that the standard deviation will
  EXCEED that threshold. Only considered if `use_admin_effect` is TRUE.

- `prior_covariate_effect`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A named list specifying the penalized complexity prior for all
  covariate effects except for the intercept, if an intercept is
  included. The two named items are "threshold", the test threshold for
  the size of each fixed effect, and "prob_above", the prior probability
  that the beta for each covariate will EXCEED that threshold. Only
  considered if `use_covariates` is TRUE and `use_stacking` is FALSE.

- `inla_link`:

  (`character(1)`)  
  Link function for fitting the INLA model, typically related to the GLM
  `family`.

- `inverse_link`:

  (`character(1)`)  
  Inverse function of `inla_link`.

- `inla_family`:

  (character)  
  GLM family to use. For more information, see
  [stats::family](https://rdrr.io/r/stats/family.html).

- `nugget_in_predict`:

  (`logical(1)`)  
  If the nugget is used in model fitting, should it also be included as
  an IID effect by pixel in the model prediction step?

- `verbose`:

  Should model progress be timed?

- `model_covariates`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  A list of covariates to be included in the INLA model. Either equal to
  `covariate_rasters`, or ML model predictions for stacked
  generalization.

- `inla_inputs_list`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  List of model inputs yielded by
  [prepare_inla_data_stack](https://henryspatialanalysis.github.io/mbg/reference/prepare_inla_data_stack.md)

- `inla_fitted_model`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  List of model outputs yielded by
  [fit_inla_model](https://henryspatialanalysis.github.io/mbg/reference/fit_inla_model.md)

- `grid_cell_predictions`:

  List of predictive surfaces yielded by
  [generate_cell_draws_and_summarize](https://henryspatialanalysis.github.io/mbg/reference/generate_cell_draws_and_summarize.md)

- `aggregated_predictions`:

  List of predictions by administrative unit. Only created if
  `aggregation_table` and `aggregation_levels` are both defined.

## Methods

### Public methods

- [`MbgModelRunner$new()`](#method-MbgModelRunner-new)

- [`MbgModelRunner$prepare_covariates()`](#method-MbgModelRunner-prepare_covariates)

- [`MbgModelRunner$fit_mbg_model()`](#method-MbgModelRunner-fit_mbg_model)

- [`MbgModelRunner$generate_predictions()`](#method-MbgModelRunner-generate_predictions)

- [`MbgModelRunner$aggregate_predictions()`](#method-MbgModelRunner-aggregate_predictions)

- [`MbgModelRunner$run_mbg_pipeline()`](#method-MbgModelRunner-run_mbg_pipeline)

- [`MbgModelRunner$get_predictive_validity()`](#method-MbgModelRunner-get_predictive_validity)

- [`MbgModelRunner$clone()`](#method-MbgModelRunner-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new MbgModelRunner object

#### Usage

    MbgModelRunner$new(
      input_data,
      id_raster,
      covariate_rasters = NULL,
      aggregation_table = NULL,
      aggregation_levels = NULL,
      population_raster = NULL,
      admin_bounds = NULL,
      admin_bounds_id = NULL,
      use_covariates = TRUE,
      use_gp = TRUE,
      use_admin_effect = FALSE,
      use_nugget = TRUE,
      use_stacking = FALSE,
      stacking_cv_settings = list(method = "repeatedcv", number = 5, repeats = 5),
      stacking_model_settings = list(gbm = NULL, treebag = NULL, rf = NULL),
      stacking_use_admin_bounds = FALSE,
      stacking_prediction_range = NULL,
      mesh_max_edge = c(0.2, 5),
      mesh_cutoff = c(0.04),
      spde_integrate_to_zero = FALSE,
      prior_spde_range = list(threshold = 0.1, prob_below = 0.05),
      prior_spde_sigma = list(threshold = 3, prob_above = 0.05),
      prior_nugget = list(threshold = 3, prob_above = 0.05),
      prior_admin_effect = list(threshold = 3, prob_above = 0.05),
      prior_covariate_effect = list(threshold = 3, prob_above = 0.05),
      inla_link = "logit",
      inverse_link = "plogis",
      inla_family = "binomial",
      nugget_in_predict = TRUE,
      verbose = TRUE
    )

#### Arguments

- `input_data`:

  ([data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html))
  Table containing at least the following fields:  

  - x (`numeric`) location x position, in the same projection as the
    `id_raster`  

  - y (`numeric`) location y position, in the same projection as the
    `id_raster`  

  - indicator (`integer`) The number of events in the population  

  - samplesize (`integer`) The total population, denominator for
    `indicator`  

- `id_raster`:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))
  raster showing the total area that will be predicted using this model

- `covariate_rasters`:

  ([`list()`](https://rdrr.io/r/base/list.html), default NULL) A list
  containing all predictor covariates. Each covariate is a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the same extent and dimensions as `id_raster`.

- `aggregation_table`:

  ([data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html))
  A table created by
  [build_aggregation_table](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md),
  linking each grid cell to one or more polygons

- `aggregation_levels`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list: for each
  named item, the name is the label for that aggregation level, and the
  value is a character vector of all fields in the original polygons to
  be used for aggregation at that level.

- `population_raster`:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))
  A raster giving population for each grid cell, to be used for
  population-weighted aggregation from grid cells to polygon boundaries.
  Should have the same dimensions as `id_raster`. If no population
  raster is passed and the results are aggregated, aggregation will be
  by simple mean rather than population-weighted mean

- `admin_bounds`:

  ([sf::sf](https://r-spatial.github.io/sf/reference/sf.html), default
  `NULL`) Polygons showing the boundaries of administrative divisions
  within the study region. Only required if `use_admin_effect` OR
  `stacking_use_admin_bounds` is `TRUE`.

- `admin_bounds_id`:

  (`character`, default `NULL`) Field containing unique identifiers for
  `admin_bounds`, if passed.

- `use_covariates`:

  (`logical(1)`, default TRUE) Should covariate effects be included in
  the predictive model?

- `use_gp`:

  (`logical(1)`, default TRUE) Should a smoothed spatial surface be
  included in the predictive model?

- `use_admin_effect`:

  (`logical(1)` default FALSE) Should IID administrative-level effects
  be included in the predictive model?

- `use_nugget`:

  (`logical(1)`, default TRUE) Should an IID effect by pixel be included
  in the predictive model?

- `use_stacking`:

  (`logical(1)`, default FALSE) Should machine learning submodels be
  trained to relate the covariate rasters with the outcome data? Only
  run if `use_covariates` is `TRUE`.

- `stacking_cv_settings`:

  ([`list()`](https://rdrr.io/r/base/list.html)) How should the stacking
  submodels be cross-validated? For more information about this term,
  see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is `TRUE`.

- `stacking_model_settings`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list of
  submodels to be run. For more information about this term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `stacking_use_admin_bounds`:

  (`logical(1)`, default FALSE) Should admin boundaries be included as
  features in the stacking submodels? For more information about this
  term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `stacking_prediction_range`:

  (`numeric(2)`, default NULL) Range of possible predictions for the
  stacking submodels. For more information about this term, see
  [run_regression_submodels](https://henryspatialanalysis.github.io/mbg/reference/run_regression_submodels.md).
  Only considered if `use_stacking` is TRUE.

- `mesh_max_edge`:

  (`numeric(2)`, default c(0.2, 5)) Maximum size of the INLA SPDE mesh
  inside (1) and outside (2) of the modeled region. Only considered if
  `use_gp` is TRUE.

- `mesh_cutoff`:

  (`numeric(1)`, default 0.04) Minimum size of the INLA mesh, usually
  reached in data-dense areas. Only considered if `use_gp` is TRUE.

- `spde_integrate_to_zero`:

  (`boolean(1)`, default FALSE) Should the 'volume' under the SPDE mesh
  integrate to zero? Only considered if `use_gp` is TRUE.

- `prior_spde_range`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list specifying
  the penalized complexity prior for the SPDE range. The two named items
  are "threshold", the test threshold (set as a proportion of the
  overall mesh extent), and "prob_below", the prior probability that the
  value is BELOW that range threshold. The function automatically
  converts "threshold" from a proportion of the overall mesh extent into
  a distance. Only considered if `use_gp` is TRUE.

- `prior_spde_sigma`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list specifying
  the penalized complexity prior for sigma (standard deviation) of the
  SPDE object. The two named items are "threshold", the test threshold
  for the standard deviation, and "prob_above", the prior probability
  that sigma will EXCEED that threshold. Only considered if `use_gp` is
  TRUE

- `prior_nugget`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list specifying
  the penalized complexity prior for the nugget term. The two named
  items are "threshold", the test threshold for the nugget standard
  deviation, and "prob_above", the prior probability that the standard
  deviation will EXCEED that threshold. Only considered if `use_nugget`
  is TRUE.

- `prior_admin_effect`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list specifying
  the penalized complexity prior for the admin-level IID term. The two
  named items are "threshold", the test threshold for the standard
  deviation of admin-level effects, and "prob_above", the prior
  probability that the standard deviation will EXCEED that threshold.
  Only considered if `use_admin_effect` is TRUE.

- `prior_covariate_effect`:

  ([`list()`](https://rdrr.io/r/base/list.html)) A named list specifying
  the penalized complexity prior for all covariate effects except for
  the intercept, if an intercept is included. The two named items are
  "threshold", the test threshold for the size of each fixed effect, and
  "prob_above", the prior probability that the beta for each covariate
  will exceed that threshold. Only considered if `use_covariates` is
  TRUE and `use_stacking` is FALSE.

- `inla_link`:

  (`character(1)`, default 'logit') Link function for fitting the INLA
  model, typically related to the GLM `family`.

- `inverse_link`:

  (`character(1)`, default 'plogis') Inverse function of `inla_link`.

- `inla_family`:

  (`character(1)`, default 'binomial') GLM family to use. For more
  information, see
  [`stats::family()`](https://rdrr.io/r/stats/family.html).

- `nugget_in_predict`:

  (`logical(1)`, default TRUE) If the nugget is used in model fitting,
  should it also be included as an IID effect by pixel in the model
  prediction step?

- `verbose`:

  (`logical(1)`, default TRUE) Should model progress be timed?

------------------------------------------------------------------------

### Method `prepare_covariates()`

Prepare covariates for MBG model fitting

#### Usage

    MbgModelRunner$prepare_covariates()

------------------------------------------------------------------------

### Method `fit_mbg_model()`

Fit MBG model

#### Usage

    MbgModelRunner$fit_mbg_model()

------------------------------------------------------------------------

### Method `generate_predictions()`

Generate predictions by grid cell

#### Usage

    MbgModelRunner$generate_predictions(n_samples = 1000, ui_width = 0.95)

#### Arguments

- `n_samples`:

  (`integer(1)`, default 1000) Number of posterior predictive samples to
  generate from the fitted model object.

- `ui_width`:

  (`numeric(1)`, default 0.95) Uncertainty interval width. This method
  will create summary rasters for quantiles ((1 - ui_width)/2) and (1 -
  (1 - ui_width)/2).

------------------------------------------------------------------------

### Method `aggregate_predictions()`

Aggregate grid cell predictions

#### Usage

    MbgModelRunner$aggregate_predictions(ui_width = 0.95)

#### Arguments

- `ui_width`:

  (`numeric(1)`, default 0.95) Uncertainty interval width. This method
  will create summary "upper" and "lower" fields rasters for quantiles
  ((1 - ui_width)/2) and (1 - (1 - ui_width)/2).

#### Returns

List with the same names as `self$aggregation_levels`, aggregating by
the columns specified in `self$aggregation_levels`

------------------------------------------------------------------------

### Method `run_mbg_pipeline()`

Run a full MBG pipeline, including stacking, MBG model fitting, and
prediction

#### Usage

    MbgModelRunner$run_mbg_pipeline(n_samples = 1000, ui_width = 0.95)

#### Arguments

- `n_samples`:

  (`integer(1)`, default 1000) Number of posterior predictive samples to
  generate from the fitted model object.

- `ui_width`:

  (`numeric(1)`, default 0.95)

------------------------------------------------------------------------

### Method `get_predictive_validity()`

Get predictive validity metrics for the fitted model

#### Usage

    MbgModelRunner$get_predictive_validity(
      in_sample = TRUE,
      validation_data = NULL,
      na.rm = FALSE
    )

#### Arguments

- `in_sample`:

  (`logical(1)`, default TRUE) Compare model predictions to the data
  used to generate the model? If FALSE, does not return the WAIC, which
  is only useful for in-sample predictive validity.

- `validation_data`:

  ([data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html),
  default NULL) Observed data to compare against. Expected for
  out-of-sample model validation. Table containing at least the
  following fields:  

  - x (`numeric`) location x position, in the same projection as the
    `id_raster`  

  - y (`numeric`) location y position, in the same projection as the
    `id_raster`  

  - indicator (`integer`) The number of events in the population  

  - samplesize (`integer`) The total population, denominator for
    `indicator`

- `na.rm`:

  (`logical(1)`, default FALSE) Should NA values be dropped from the
  RMSE and log predictive density calculations?

#### Details

Returns the point RMSE (compared against the mean estimates by pixel),
log-posterior density (compared against the predictive draws), and the
Watanabe-Aikake Information Criterion (WAIC, only returned for in-sample
predictive validity).

#### Returns

[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
Containing the following fields:  

- 'rmse': Root mean squared error when compared against the mean
  estimates by pixel. Lower RMSE is better.  

- 'lpd': Log posterior predictive density when compared against
  pixel-level samples from the model. Higher LPD is better.  

- 'waic' (in-sample only): Watanable-Aikake information criterion
  estimated by INLA. Lower WAIC is better.  
  For clarity, these fields will have the suffix "\_is" for in-sample
  models, and "\_oos" for out-of-sample models.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MbgModelRunner$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
