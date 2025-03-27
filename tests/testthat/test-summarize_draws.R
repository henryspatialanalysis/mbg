# Create test data
test_matrix <- matrix(
  c(1, 2, 3, 4, 5,
    2, 3, 4, 5, 6,
    3, 4, 5, 6, 7),
  nrow = 3,
  ncol = 5
)

test_df <- data.frame(
  id = c("A", "B", "C"),
  pop = c(100, 200, 300),
  draw1 = c(1, 2, 3),
  draw2 = c(2, 3, 4),
  draw3 = c(3, 4, 5)
)

test_that("summarize_draws works with matrix input", {
  # Test basic functionality
  result <- summarize_draws(test_matrix)
  expect_true(inherits(result, "data.table"))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)  # mean, lower, upper, ui_width

  # Test with different UI width
  result_wide <- summarize_draws(test_matrix, ui_width = 0.99)
  expect_true(result_wide$ui_width > result$ui_width)

  # Test with na.rm = FALSE
  na_matrix <- test_matrix
  na_matrix[1, 1] <- NA
  result_na <- summarize_draws(na_matrix, na.rm = FALSE)
  expect_true(is.na(result_na$mean[1]))
})

test_that("summarize_draws works with data.frame input", {
  # Test basic functionality
  result <- summarize_draws(test_df)
  expect_true(inherits(result, "data.table"))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)  # mean, lower, upper, ui_width

  # Test with id_fields
  result_id <- summarize_draws(test_df, id_fields = "id")
  expect_true(inherits(result_id, "data.table"))
  expect_equal(nrow(result_id), 3)
  expect_equal(ncol(result_id), 5)  # id, mean, lower, upper, ui_width
  expect_true("id" %in% colnames(result_id))

  # Test with draw_fields
  result_draw <- summarize_draws(test_df, draw_fields = c("draw1", "draw2"))
  expect_true(inherits(result_draw, "data.table"))
  expect_equal(nrow(result_draw), 3)
  expect_equal(ncol(result_draw), 4)  # mean, lower, upper, ui_width
})

test_that("summarize_draws handles edge cases", {
  # Test with single draw
  single_draw <- data.frame(
    id = c("A", "B"),
    draw1 = c(1, 2)
  )
  result_single <- summarize_draws(single_draw, id_fields = "id")
  expect_true(inherits(result_single, "data.table"))
  expect_equal(nrow(result_single), 2)
  expect_equal(result_single$mean, result_single$lower)
  expect_equal(result_single$mean, result_single$upper)

  # Test with all NA values
  na_df <- data.frame(
    id = c("A", "B"),
    draw1 = c(NA, NA),
    draw2 = c(NA, NA)
  )
  result_na <- summarize_draws(na_df, id_fields = "id")
  expect_true(all(is.na(result_na$mean)))
  expect_true(all(is.na(result_na$lower)))
  expect_true(all(is.na(result_na$upper)))
  expect_true(all(is.na(result_na$ui_width)))
})

test_that("summarize_draws validates inputs correctly", {
  # Test invalid draws type
  expect_error(summarize_draws("not a matrix or data frame"))

  # Test invalid id_fields
  expect_error(summarize_draws(test_df, id_fields = "nonexistent"))

  # Test invalid draw_fields
  expect_error(summarize_draws(test_df, draw_fields = "nonexistent"))

  # Test invalid ui_width
  expect_error(summarize_draws(test_matrix, ui_width = "not numeric"))
  expect_error(summarize_draws(test_matrix, ui_width = 0))
  expect_error(summarize_draws(test_matrix, ui_width = 1.5))

  # Test invalid na.rm
  expect_error(summarize_draws(test_matrix, na.rm = "not logical"))
})

test_that("summarize_draws preserves data types", {
  # Test with integer draws
  int_df <- data.frame(
    id = c("A", "B"),
    draw1 = c(1L, 2L),
    draw2 = c(2L, 3L)
  )
  result_int <- summarize_draws(int_df, id_fields = "id")
  expect_true(is.character(result_int$id))
  expect_true(is.numeric(result_int$mean))

  # Test with numeric draws
  num_df <- data.frame(
    id = c("A", "B"),
    draw1 = c(1.5, 2.5),
    draw2 = c(2.5, 3.5)
  )
  result_num <- summarize_draws(num_df, id_fields = "id")
  expect_true(is.character(result_num$id))
  expect_true(is.numeric(result_num$mean))
})
