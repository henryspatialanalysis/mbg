# Create test data
test_data <- data.frame(
  x1 = stats::rnorm(100),
  x2 = stats::rnorm(100),
  x3 = stats::rnorm(100),
  x4 = stats::rnorm(100)
)

# Create correlated data
correlated_data <- data.frame(
  x1 = stats::rnorm(100),
  x2 = x1 + stats::rnorm(100, sd = 0.1),  # Highly correlated with x1
  x3 = stats::rnorm(100),
  x4 = x3 + stats::rnorm(100, sd = 0.1)   # Highly correlated with x3
)

test_that("vif_covariate_select works correctly", {
  # Test with uncorrelated data
  result <- vif_covariate_select(test_data)
  expect_true(inherits(result, "data.table"))
  expect_equal(nrow(result), 4)
  expect_true(all(result$keep))  # All variables should be kept

  # Test with correlated data
  result_corr <- vif_covariate_select(correlated_data)
  expect_true(inherits(result_corr, "data.table"))
  expect_equal(nrow(result_corr), 4)
  expect_false(all(result_corr$keep))  # Some variables should be dropped

  # Test with different cutoff
  result_strict <- vif_covariate_select(correlated_data, vif_cutoff = 2)
  expect_true(inherits(result_strict, "data.table"))
  expect_equal(nrow(result_strict), 4)
  expect_false(all(result_strict$keep))  # More variables should be dropped
})

test_that("vif_covariate_select handles edge cases", {
  # Test with single column
  single_col <- data.frame(x = stats::rnorm(100))
  result_single <- vif_covariate_select(single_col)
  expect_true(inherits(result_single, "data.table"))
  expect_equal(nrow(result_single), 1)
  expect_true(result_single$keep)

  # Test with two columns
  two_cols <- data.frame(x1 = stats::rnorm(100), x2 = stats::rnorm(100))
  result_two <- vif_covariate_select(two_cols)
  expect_true(inherits(result_two, "data.table"))
  expect_equal(nrow(result_two), 2)
  expect_true(all(result_two$keep))

  # Test with three columns
  three_cols <- data.frame(
    x1 = stats::rnorm(100),
    x2 = x1 + stats::rnorm(100, sd = 0.1),
    x3 = stats::rnorm(100)
  )
  result_three <- vif_covariate_select(three_cols)
  expect_true(inherits(result_three, "data.table"))
  expect_equal(nrow(result_three), 3)
  expect_false(all(result_three$keep))
})

test_that("vif_covariate_select validates inputs correctly", {
  # Test invalid dataset type
  expect_error(vif_covariate_select("not a data frame"))

  # Test duplicate column names
  dup_cols <- data.frame(
    x = stats::rnorm(100),
    x = stats::rnorm(100)  # Duplicate column name
  )
  expect_error(vif_covariate_select(dup_cols))

  # Test invalid vif_cutoff
  expect_error(vif_covariate_select(test_data, vif_cutoff = "not numeric"))
  expect_error(vif_covariate_select(test_data, vif_cutoff = NA))
  expect_error(vif_covariate_select(test_data, vif_cutoff = 0))
  expect_error(vif_covariate_select(test_data, vif_cutoff = -1))
})

test_that("vif_covariate_select handles missing values", {
  # Create data with missing values
  na_data <- test_data
  na_data$x1[1:10] <- NA

  # Test with na.rm = TRUE (default)
  result_na <- vif_covariate_select(na_data)
  expect_true(inherits(result_na, "data.table"))
  expect_equal(nrow(result_na), 4)

  # Test with complete data
  complete_data <- stats::na.omit(na_data)
  result_complete <- vif_covariate_select(complete_data)
  expect_true(inherits(result_complete, "data.table"))
  expect_equal(nrow(result_complete), 4)
})