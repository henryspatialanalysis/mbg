# Create test data
test_poly <- terra::vect("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))")
test_id_rast <- build_id_raster(test_poly)

# Create test values
test_values <- 1:4  # 4 non-NA pixels in test_id_rast
test_matrix <- matrix(1:8, nrow = 4, ncol = 2)  # 2 layers
test_df <- data.frame(
  layer1 = 1:4,
  layer2 = 5:8
)

test_that("values_to_raster works with vector input", {
  # Test basic functionality
  result <- values_to_raster(test_values, test_id_rast)
  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 1)
  expect_equal(terra::values(result)[!is.na(terra::values(result))], test_values)

  # Test with multiple layers
  multi_values <- rep(test_values, 2)  # Two layers
  result_multi <- values_to_raster(multi_values, test_id_rast)
  expect_true(inherits(result_multi, "SpatRaster"))
  expect_equal(terra::nlyr(result_multi), 2)
})

test_that("values_to_raster works with matrix input", {
  # Test basic functionality
  result <- values_to_raster(test_matrix, test_id_rast)
  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 2)
  expect_equal(terra::values(result)[!is.na(terra::values(result))], as.vector(test_matrix))
})

test_that("values_to_raster works with data.frame input", {
  # Test basic functionality
  result <- values_to_raster(test_df, test_id_rast)
  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 2)
  expect_equal(terra::values(result)[!is.na(terra::values(result))], as.vector(as.matrix(test_df)))
})

test_that("values_to_raster validates inputs correctly", {
  # Test invalid x type
  expect_error(values_to_raster("not numeric", test_id_rast))

  # Test invalid id_raster type
  expect_error(values_to_raster(test_values, "not a raster"))

  # Test id_raster with multiple layers
  multi_layer_rast <- c(test_id_rast, test_id_rast)
  expect_error(values_to_raster(test_values, multi_layer_rast))

  # Test length mismatch
  expect_error(values_to_raster(1:5, test_id_rast))  # Wrong number of values
})

test_that("values_to_raster handles edge cases", {
  # Test with empty vector
  expect_error(values_to_raster(numeric(0), test_id_rast))

  # Test with all NA values
  na_rast <- test_id_rast
  terra::values(na_rast) <- NA
  expect_error(values_to_raster(test_values, na_rast))

  # Test with single value
  single_value <- 1
  result_single <- values_to_raster(single_value, test_id_rast)
  expect_true(inherits(result_single, "SpatRaster"))
  expect_equal(terra::nlyr(result_single), 1)
  expect_true(all(terra::values(result_single)[!is.na(terra::values(result_single))] == 1))
})

test_that("values_to_raster preserves data types", {
  # Test with integer values
  int_values <- 1L:4L
  result_int <- values_to_raster(int_values, test_id_rast)
  expect_true(is.integer(terra::values(result_int)[!is.na(terra::values(result_int))]))

  # Test with numeric values
  num_values <- c(1.1, 2.2, 3.3, 4.4)
  result_num <- values_to_raster(num_values, test_id_rast)
  expect_true(is.numeric(terra::values(result_num)[!is.na(terra::values(result_num))]))
})
