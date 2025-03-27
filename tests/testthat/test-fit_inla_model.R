# Create test data
test_poly <- terra::vect("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))")
test_id_rast <- build_id_raster(test_poly)

# Create test input data
test_input_data <- data.table::data.table(
  x = c(0.5, 1.5),
  y = c(0.5, 1.5),
  indicator = c(1, 0),
  samplesize = c(10, 10),
  cluster_id = c(1, 2)
)

# Create test covariates
test_covariate <- terra::rast(matrix(1:4, 2, 2))
test_covariates <- list(cov1 = test_covariate)

# Prepare INLA data stack
test_stack <- prepare_inla_data_stack(
  input_data = test_input_data,
  id_raster = test_id_rast,
  covariates = test_covariates
)

test_that("fit_inla_model works with basic inputs", {
  # Test basic functionality
  result <- fit_inla_model(
    formula = test_stack$formula_string,
    data_stack = test_stack$inla_data_stack,
    spde = test_stack$spde
  )
  expect_true(inherits(result, "inla"))

  # Test with different family
  result_gaussian <- fit_inla_model(
    formula = test_stack$formula_string,
    data_stack = test_stack$inla_data_stack,
    spde = test_stack$spde,
    family = "gaussian",
    link = "identity"
  )
  expect_true(inherits(result_gaussian, "inla"))
})

test_that("fit_inla_model handles samplesize and precision correctly", {
  # Test with custom samplesize
  result_samplesize <- fit_inla_model(
    formula = test_stack$formula_string,
    data_stack = test_stack$inla_data_stack,
    spde = test_stack$spde,
    samplesize_vec = c(20, 20)
  )
  expect_true(inherits(result_samplesize, "inla"))

  # Test with custom precision
  result_precision <- fit_inla_model(
    formula = test_stack$formula_string,
    data_stack = test_stack$inla_data_stack,
    spde = test_stack$spde,
    family = "gaussian",
    link = "identity",
    precision_vec = c(0.5, 0.5)
  )
  expect_true(inherits(result_precision, "inla"))
})

test_that("fit_inla_model handles fixed effects priors correctly", {
  # Test with custom fixed effects prior
  result_prior <- fit_inla_model(
    formula = test_stack$formula_string,
    data_stack = test_stack$inla_data_stack,
    spde = test_stack$spde,
    fixed_effects_pc_prior = list(threshold = 5, prob_above = 0.1)
  )
  expect_true(inherits(result_prior, "inla"))
})

test_that("fit_inla_model handles verbose output correctly", {
  # Test with verbose = TRUE
  expect_silent(
    result_verbose <- fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      verbose = TRUE
    )
  )
  expect_true(inherits(result_verbose, "inla"))

  # Test with verbose = FALSE
  expect_silent(
    result_quiet <- fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      verbose = FALSE
    )
  )
  expect_true(inherits(result_quiet, "inla"))
})

test_that("fit_inla_model validates inputs correctly", {
  # Test invalid formula
  expect_error(
    fit_inla_model(
      formula = "not a formula",
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde
    )
  )

  # Test invalid data stack
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = "not a stack",
      spde = test_stack$spde
    )
  )

  # Test invalid spde
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = "not an spde"
    )
  )

  # Test invalid samplesize
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      samplesize_vec = c(1, 2, 3)  # Wrong length
    )
  )

  # Test invalid precision
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      family = "gaussian",
      precision_vec = c(1, 2, 3)  # Wrong length
    )
  )

  # Test invalid family
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      family = "invalid"
    )
  )

  # Test invalid link
  expect_error(
    fit_inla_model(
      formula = test_stack$formula_string,
      data_stack = test_stack$inla_data_stack,
      spde = test_stack$spde,
      link = "invalid"
    )
  )
})

test_that("fit_inla_model handles edge cases", {
  # Test with single observation
  single_data <- test_input_data[1]
  single_stack <- prepare_inla_data_stack(
    input_data = single_data,
    id_raster = test_id_rast,
    covariates = test_covariates
  )
  result_single <- fit_inla_model(
    formula = single_stack$formula_string,
    data_stack = single_stack$inla_data_stack,
    spde = single_stack$spde
  )
  expect_true(inherits(result_single, "inla"))

  # Test with no covariates
  no_cov_stack <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_covariates = FALSE
  )
  result_no_cov <- fit_inla_model(
    formula = no_cov_stack$formula_string,
    data_stack = no_cov_stack$inla_data_stack,
    spde = no_cov_stack$spde
  )
  expect_true(inherits(result_no_cov, "inla"))
})
