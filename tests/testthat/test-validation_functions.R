# Create test data
test_rast <- terra::rast(matrix(1:4, 2, 2))
test_id_rast <- terra::rast(matrix(1:4, 2, 2))
test_validation_data <- data.frame(
  x = c(0.5, 1.5),
  y = c(0.5, 1.5),
  indicator = c(1, 0),
  samplesize = c(10, 10)
)

test_that("rmse_raster_to_point works correctly", {
  # Test basic functionality
  result <- rmse_raster_to_point(
    estimates = test_rast,
    validation_data = test_validation_data,
    outcome_field = "indicator"
  )
  expect_true(is.numeric(result))
  expect_true(length(result) == 1)

  # Test with na.rm = TRUE
  na_data <- test_validation_data
  na_data$indicator[1] <- NA
  result_na <- rmse_raster_to_point(
    estimates = test_rast,
    validation_data = na_data,
    outcome_field = "indicator",
    na.rm = TRUE
  )
  expect_true(is.numeric(result_na))

  # Test with na.rm = FALSE
  result_na_false <- rmse_raster_to_point(
    estimates = test_rast,
    validation_data = na_data,
    outcome_field = "indicator",
    na.rm = FALSE
  )
  expect_true(is.na(result_na_false))
})

test_that("rmse_raster_to_point validates inputs correctly", {
  # Test invalid outcome field
  expect_error(
    rmse_raster_to_point(
      estimates = test_rast,
      validation_data = test_validation_data,
      outcome_field = "nonexistent"
    )
  )

  # Test invalid coordinates
  invalid_data <- test_validation_data
  invalid_data$x[1] <- NA
  expect_error(
    rmse_raster_to_point(
      estimates = test_rast,
      validation_data = invalid_data,
      outcome_field = "indicator"
    )
  )

  # Test invalid na.rm parameter
  expect_error(
    rmse_raster_to_point(
      estimates = test_rast,
      validation_data = test_validation_data,
      outcome_field = "indicator",
      na.rm = "not logical"
    )
  )
})

test_that("log_posterior_density works correctly", {
  # Create test draws matrix
  test_draws <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 4,
    ncol = 1
  )

  # Test basic functionality
  result <- log_posterior_density(
    draws = test_draws,
    validation_data = test_validation_data,
    id_raster = test_id_rast
  )
  expect_true(is.numeric(result))
  expect_true(length(result) == 1)

  # Test with na.rm = TRUE
  na_data <- test_validation_data
  na_data$indicator[1] <- NA
  result_na <- log_posterior_density(
    draws = test_draws,
    validation_data = na_data,
    id_raster = test_id_rast,
    na.rm = TRUE
  )
  expect_true(is.numeric(result_na))

  # Test with na.rm = FALSE
  result_na_false <- log_posterior_density(
    draws = test_draws,
    validation_data = na_data,
    id_raster = test_id_rast,
    na.rm = FALSE
  )
  expect_true(is.na(result_na_false))
})

test_that("log_posterior_density validates inputs correctly", {
  # Create test draws matrix
  test_draws <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 4,
    ncol = 1
  )

  # Test invalid draws matrix
  invalid_draws <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.5),
    nrow = 5,
    ncol = 1
  )
  expect_error(
    log_posterior_density(
      draws = invalid_draws,
      validation_data = test_validation_data,
      id_raster = test_id_rast
    )
  )

  # Test invalid validation data
  invalid_data <- test_validation_data
  invalid_data$indicator[1] <- -1  # Invalid binomial count
  expect_error(
    log_posterior_density(
      draws = test_draws,
      validation_data = invalid_data,
      id_raster = test_id_rast
    )
  )

  # Test invalid coordinates
  invalid_coords <- test_validation_data
  invalid_coords$x[1] <- NA
  expect_error(
    log_posterior_density(
      draws = test_draws,
      validation_data = invalid_coords,
      id_raster = test_id_rast
    )
  )

  # Test invalid na.rm parameter
  expect_error(
    log_posterior_density(
      draws = test_draws,
      validation_data = test_validation_data,
      id_raster = test_id_rast,
      na.rm = "not logical"
    )
  )
})