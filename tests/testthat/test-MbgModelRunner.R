# Create test data
test_rast <- terra::rast(matrix(1:4, 2, 2))
test_poly <- terra::vect("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))")
test_agg_table <- build_aggregation_table(
  data_raster = test_rast,
  polygon_sf = test_poly,
  id_col = "polygon_id"
)

# Create test input data
test_input_data <- data.table::data.table(
  x = c(0.5, 1.5),
  y = c(0.5, 1.5),
  indicator = c(1, 0),
  samplesize = c(10, 10)
)

# Create test covariate raster
test_covariate <- terra::rast(matrix(c(1, 2, 3, 4), 2, 2))
test_covariates <- list(cov1 = test_covariate)

test_that("MbgModelRunner initializes correctly", {
  # Test basic initialization
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast
  )
  expect_true(inherits(runner, "MbgModelRunner"))
  expect_equal(runner$input_data, test_input_data)
  expect_equal(runner$id_raster, test_rast)

  # Test initialization with covariates
  runner_cov <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates
  )
  expect_equal(runner_cov$covariate_rasters, test_covariates)

  # Test initialization with aggregation
  runner_agg <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    aggregation_table = test_agg_table,
    aggregation_levels = list(admin0 = "polygon_id")
  )
  expect_equal(runner_agg$aggregation_table, test_agg_table)
  expect_equal(runner_agg$aggregation_levels, list(admin0 = "polygon_id"))
})

test_that("MbgModelRunner validates inputs correctly", {
  # Test invalid input data
  invalid_data <- data.table::data.table(
    x = c(0.5, 1.5),
    y = c(0.5, 1.5),
    indicator = c(1, 0)
    # Missing samplesize
  )
  expect_error(
    MbgModelRunner$new(
      input_data = invalid_data,
      id_raster = test_rast
    )
  )

  # Test invalid raster dimensions
  invalid_rast <- terra::rast(matrix(1:6, 2, 3))  # Different dimensions
  expect_error(
    MbgModelRunner$new(
      input_data = test_input_data,
      id_raster = test_rast,
      covariate_rasters = list(cov1 = invalid_rast)
    )
  )

  # Test invalid aggregation table
  invalid_table <- data.table::data.table(pixel_id = 1:5)  # Missing required columns
  expect_error(
    MbgModelRunner$new(
      input_data = test_input_data,
      id_raster = test_rast,
      aggregation_table = invalid_table
    )
  )
})

test_that("MbgModelRunner handles model fitting correctly", {
  # Create a basic model runner
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates,
    use_covariates = TRUE,
    use_gp = TRUE,
    use_nugget = TRUE,
    inla_family = "binomial",
    inla_link = "logit",
    inverse_link = "inv.logit"
  )

  # Test model fitting
  expect_silent(runner$fit_model())
  expect_true(!is.null(runner$inla_fitted_model))
  expect_true(!is.null(runner$inla_inputs_list))

  # Test model fitting with different settings
  runner_no_gp <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates,
    use_covariates = TRUE,
    use_gp = FALSE,
    use_nugget = TRUE,
    inla_family = "binomial",
    inla_link = "logit",
    inverse_link = "inv.logit"
  )
  expect_silent(runner_no_gp$fit_model())
})

test_that("MbgModelRunner handles predictions correctly", {
  # Create and fit a model
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates,
    use_covariates = TRUE,
    use_gp = TRUE,
    use_nugget = TRUE,
    inla_family = "binomial",
    inla_link = "logit",
    inverse_link = "inv.logit"
  )
  runner$fit_model()

  # Test grid cell predictions
  expect_silent(runner$predict_grid_cells())
  expect_true(!is.null(runner$grid_cell_predictions))

  # Test aggregation
  runner_agg <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates,
    aggregation_table = test_agg_table,
    aggregation_levels = list(admin0 = "polygon_id"),
    use_covariates = TRUE,
    use_gp = TRUE,
    use_nugget = TRUE,
    inla_family = "binomial",
    inla_link = "logit",
    inverse_link = "inv.logit"
  )
  runner_agg$fit_model()
  runner_agg$predict_grid_cells()
  expect_silent(runner_agg$aggregate_predictions())
  expect_true(!is.null(runner_agg$aggregated_predictions))
})

test_that("MbgModelRunner handles stacking correctly", {
  # Create a model runner with stacking
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    covariate_rasters = test_covariates,
    use_covariates = TRUE,
    use_stacking = TRUE,
    stacking_model_settings = list(
      gbm = list(
        method = "gbm",
        tuneLength = 2
      )
    ),
    stacking_cv_settings = list(
      V = 2,
      stratify = FALSE
    ),
    inla_family = "binomial",
    inla_link = "logit",
    inverse_link = "inv.logit"
  )

  # Test stacking model fitting
  expect_silent(runner$fit_stacking_models())
  expect_true(!is.null(runner$model_covariates))

  # Test full model fitting with stacking
  expect_silent(runner$fit_model())
  expect_true(!is.null(runner$inla_fitted_model))
})

test_that("MbgModelRunner handles errors gracefully", {
  # Test prediction before fitting
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast
  )
  expect_error(runner$predict_grid_cells())

  # Test aggregation without aggregation table
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast
  )
  runner$fit_model()
  runner$predict_grid_cells()
  expect_error(runner$aggregate_predictions())

  # Test stacking without covariates
  runner <- MbgModelRunner$new(
    input_data = test_input_data,
    id_raster = test_rast,
    use_covariates = FALSE,
    use_stacking = TRUE
  )
  expect_error(runner$fit_stacking_models())
})
