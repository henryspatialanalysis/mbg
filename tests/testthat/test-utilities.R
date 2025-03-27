# Create test sf object
test_poly <- sf::st_as_sf(data.frame(
  id = c(1, 2, 3),
  group = c("A", "A", "B")
), geometry = sf::st_sfc(
  sf::st_polygon(list(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0)))),
  sf::st_polygon(list(rbind(c(1,0), c(1,1), c(2,1), c(2,0), c(1,0)))),
  sf::st_polygon(list(rbind(c(2,0), c(2,1), c(3,1), c(3,0), c(2,0))))
))

test_that("make_time_stamp works correctly", {
  # Test basic functionality
  result <- make_time_stamp()
  expect_true(is.character(result))
  expect_true(nchar(result) >= 15)  # YYYYMMDD_HH_MM_SS format

  # Test with suffix
  result_suffix <- make_time_stamp(suffix = "test")
  expect_true(endsWith(result_suffix, "_test"))

  # Test without milliseconds
  result_no_ms <- make_time_stamp(milliseconds = FALSE)
  expect_true(nchar(result_no_ms) == 15)  # YYYYMMDD_HH_MM_SS format

  # Test with both suffix and no milliseconds
  result_combo <- make_time_stamp(suffix = "test", milliseconds = FALSE)
  expect_true(endsWith(result_combo, "_test"))
  expect_true(nchar(result_combo) == 20)  # YYYYMMDD_HH_MM_SS_test format
})

test_that("make_time_stamp validates inputs correctly", {
  # Test invalid suffix length
  expect_error(make_time_stamp(suffix = c("test1", "test2")))

  # Test invalid milliseconds parameter
  expect_error(make_time_stamp(milliseconds = "not logical"))
})

test_that("dissolve_sf_by_attribute works correctly", {
  # Test dissolving all polygons
  result_all <- dissolve_sf_by_attribute(test_poly)
  expect_true(inherits(result_all, "sf"))
  expect_equal(nrow(result_all), 1)

  # Test dissolving by attribute
  result_by_group <- dissolve_sf_by_attribute(test_poly, by = "group")
  expect_true(inherits(result_by_group, "sf"))
  expect_equal(nrow(result_by_group), 2)  # Two unique groups

  # Test dissolving by multiple attributes
  test_poly$subgroup <- c("X", "X", "Y")
  result_by_multiple <- dissolve_sf_by_attribute(
    test_poly,
    by = c("group", "subgroup")
  )
  expect_true(inherits(result_by_multiple, "sf"))
  expect_equal(nrow(result_by_multiple), 3)  # Three unique combinations
})

test_that("dissolve_sf_by_attribute validates inputs correctly", {
  # Create test sf object
  test_poly <- sf::st_as_sf(data.frame(
    id = c(1, 2),
    group = c("A", "B")
  ), geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0)))),
    sf::st_polygon(list(rbind(c(1,0), c(1,1), c(2,1), c(2,0), c(1,0))))
  ))

  # Test invalid input type
  expect_error(dissolve_sf_by_attribute("not an sf object"))

  # Test invalid by parameter
  expect_error(dissolve_sf_by_attribute(test_poly, by = 1))  # Not character
  expect_error(dissolve_sf_by_attribute(test_poly, by = "nonexistent"))  # Column doesn't exist
})