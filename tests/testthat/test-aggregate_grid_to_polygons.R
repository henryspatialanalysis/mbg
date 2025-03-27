# Create test data
test_rast <- terra::rast(matrix(1:4, 2, 2))
test_poly <- terra::vect("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))")
test_agg_table <- build_aggregation_table(
  data_raster = test_rast,
  polygon_sf = test_poly,
  id_col = "polygon_id"
)

test_that("aggregate_raster_to_polygons works with basic inputs", {
  # Test mean aggregation
  result <- aggregate_raster_to_polygons(
    data_raster = test_rast,
    aggregation_table = test_agg_table,
    method = "mean"
  )
  expect_true(inherits(result, "data.table"))
  expect_equal(nrow(result), 1)
  expect_equal(result$data, 2.5)

  # Test sum aggregation
  result_sum <- aggregate_raster_to_polygons(
    data_raster = test_rast,
    aggregation_table = test_agg_table,
    method = "sum"
  )
  expect_equal(result_sum$data, 10)
})

test_that("aggregate_raster_to_polygons handles multiple layers", {
  # Create multi-layer raster
  multi_rast <- c(test_rast, test_rast * 2)

  # Test without z_dimension
  result <- aggregate_raster_to_polygons(
    data_raster = multi_rast,
    aggregation_table = test_agg_table
  )
  expect_equal(nrow(result), 2)
  expect_equal(result$z, c(1, 2))

  # Test with custom z_dimension
  result_custom <- aggregate_raster_to_polygons(
    data_raster = multi_rast,
    aggregation_table = test_agg_table,
    z_dimension = c("year1", "year2"),
    z_dimension_name = "year"
  )
  expect_equal(result_custom$year, c("year1", "year2"))
})

test_that("aggregate_raster_to_polygons handles weighted aggregation", {
  # Create weighting raster
  weight_rast <- terra::rast(matrix(c(1, 2, 3, 4), 2, 2))

  # Test weighted mean
  result <- aggregate_raster_to_polygons(
    data_raster = test_rast,
    aggregation_table = test_agg_table,
    method = "weighted.mean",
    weighting_raster = weight_rast
  )
  expect_true(inherits(result, "data.table"))

  # Test weighted sum
  result_sum <- aggregate_raster_to_polygons(
    data_raster = test_rast,
    aggregation_table = test_agg_table,
    method = "weighted.sum",
    weighting_raster = weight_rast
  )
  expect_true(inherits(result_sum, "data.table"))
})

test_that("aggregate_raster_to_polygons handles NA values", {
  # Create raster with NAs
  na_rast <- test_rast
  na_rast[1] <- NA

  # Test with na.rm = TRUE
  result <- aggregate_raster_to_polygons(
    data_raster = na_rast,
    aggregation_table = test_agg_table,
    na.rm = TRUE
  )
  expect_true(!is.na(result$data))

  # Test with na.rm = FALSE
  result_na <- aggregate_raster_to_polygons(
    data_raster = na_rast,
    aggregation_table = test_agg_table,
    na.rm = FALSE
  )
  expect_true(is.na(result_na$data))
})

test_that("aggregate_raster_to_polygons validates inputs correctly", {
  # Test invalid method
  expect_error(
    aggregate_raster_to_polygons(
      data_raster = test_rast,
      aggregation_table = test_agg_table,
      method = "invalid"
    )
  )

  # Test missing weighting raster for weighted methods
  expect_error(
    aggregate_raster_to_polygons(
      data_raster = test_rast,
      aggregation_table = test_agg_table,
      method = "weighted.mean"
    )
  )

  # Test invalid aggregation table
  invalid_table <- data.table::data.table(pixel_id = 1:5)  # Missing required columns
  expect_error(
    aggregate_raster_to_polygons(
      data_raster = test_rast,
      aggregation_table = invalid_table
    )
  )

  # Test invalid z_dimension length
  expect_error(
    aggregate_raster_to_polygons(
      data_raster = test_rast,
      aggregation_table = test_agg_table,
      z_dimension = 1:3  # Wrong length
    )
  )
})

test_that("aggregate_draws_to_polygons works with matrix inputs", {
  # Create test draws matrix
  draws_matrix <- matrix(1:4, 2, 2)

  # Test basic aggregation
  result <- aggregate_draws_to_polygons(
    draws_matrix = draws_matrix,
    aggregation_table = test_agg_table
  )
  expect_true(inherits(result, "data.table"))
  expect_equal(nrow(result), 1)

  # Test with multiple draws
  multi_draws <- array(1:8, dim = c(2, 2, 2))
  result_multi <- aggregate_draws_to_polygons(
    draws_matrix = multi_draws,
    aggregation_table = test_agg_table
  )
  expect_equal(nrow(result_multi), 2)
})

test_that("aggregate_draws_to_polygons handles weighted aggregation", {
  # Create test data
  draws_matrix <- matrix(1:4, 2, 2)
  weight_rast <- terra::rast(matrix(c(1, 2, 3, 4), 2, 2))

  # Test weighted mean
  result <- aggregate_draws_to_polygons(
    draws_matrix = draws_matrix,
    aggregation_table = test_agg_table,
    method = "weighted.mean",
    weighting_raster = weight_rast
  )
  expect_true(inherits(result, "data.table"))

  # Test weighted sum
  result_sum <- aggregate_draws_to_polygons(
    draws_matrix = draws_matrix,
    aggregation_table = test_agg_table,
    method = "weighted.sum",
    weighting_raster = weight_rast
  )
  expect_true(inherits(result_sum, "data.table"))
})

test_that("aggregate_draws_to_polygons validates inputs correctly", {
  # Test invalid method
  expect_error(
    aggregate_draws_to_polygons(
      draws_matrix = matrix(1:4, 2, 2),
      aggregation_table = test_agg_table,
      method = "invalid"
    )
  )

  # Test missing weighting raster for weighted methods
  expect_error(
    aggregate_draws_to_polygons(
      draws_matrix = matrix(1:4, 2, 2),
      aggregation_table = test_agg_table,
      method = "weighted.mean"
    )
  )

  # Test invalid draws matrix dimensions
  expect_error(
    aggregate_draws_to_polygons(
      draws_matrix = matrix(1:6, 2, 3),  # Wrong dimensions
      aggregation_table = test_agg_table
    )
  )
})