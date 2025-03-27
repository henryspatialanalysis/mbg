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

# Create test admin boundaries
test_admin <- sf::st_as_sf(data.frame(
  admin_id = 1:2
), geometry = sf::st_sfc(
  sf::st_polygon(list(rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0)))),
  sf::st_polygon(list(rbind(c(1,0), c(1,1), c(2,1), c(2,0), c(1,0))))
))

test_that("prepare_inla_data_stack works with basic inputs", {
  # Test basic functionality
  result <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates
  )
  expect_true(is.list(result))
  expect_true(all(c("mesh", "spde", "inla_data_stack", "formula_string") %in% names(result)))

  # Test without covariates
  result_no_cov <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_covariates = FALSE
  )
  expect_true(is.list(result_no_cov))
  expect_true(all(c("mesh", "spde", "inla_data_stack", "formula_string") %in% names(result_no_cov)))
})

test_that("prepare_inla_data_stack handles SPDE correctly", {
  # Test with SPDE
  result_spde <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_spde = TRUE
  )
  expect_true(!is.null(result_spde$mesh))
  expect_true(!is.null(result_spde$spde))

  # Test without SPDE
  result_no_spde <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_spde = FALSE
  )
  expect_true(is.null(result_no_spde$mesh))
  expect_true(is.null(result_no_spde$spde))
})

test_that("prepare_inla_data_stack handles nugget correctly", {
  # Test with nugget
  result_nugget <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_nugget = TRUE
  )
  expect_true(grepl("nugget", result_nugget$formula_string))

  # Test without nugget
  result_no_nugget <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_nugget = FALSE
  )
  expect_false(grepl("nugget", result_no_nugget$formula_string))
})

test_that("prepare_inla_data_stack handles admin effects correctly", {
  # Test with admin effects
  result_admin <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_admin_effect = TRUE,
    admin_boundaries = test_admin
  )
  expect_true(grepl("adm_effect", result_admin$formula_string))

  # Test without admin effects
  result_no_admin <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    use_admin_effect = FALSE
  )
  expect_false(grepl("adm_effect", result_no_admin$formula_string))
})

test_that("prepare_inla_data_stack handles covariates correctly", {
  # Test with covariates sum to one
  result_sum <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates,
    covariates_sum_to_one = TRUE
  )
  expect_true(grepl("extraconstr", result_sum$formula_string))

  # Test with multiple covariates
  test_covariates_multi <- list(
    cov1 = test_covariate,
    cov2 = test_covariate * 2
  )
  result_multi <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = test_covariates_multi
  )
  expect_true(grepl("covariates", result_multi$formula_string))
})

test_that("prepare_inla_data_stack validates inputs correctly", {
  # Test invalid input data
  invalid_data <- test_input_data
  invalid_data$x[1] <- NA
  expect_error(
    prepare_inla_data_stack(
      input_data = invalid_data,
      id_raster = test_id_rast,
      covariates = test_covariates
    )
  )

  # Test invalid covariates
  invalid_cov <- list(cov1 = "not a raster")
  expect_error(
    prepare_inla_data_stack(
      input_data = test_input_data,
      id_raster = test_id_rast,
      covariates = invalid_cov
    )
  )

  # Test invalid admin boundaries
  expect_error(
    prepare_inla_data_stack(
      input_data = test_input_data,
      id_raster = test_id_rast,
      covariates = test_covariates,
      use_admin_effect = TRUE,
      admin_boundaries = "not an sf object"
    )
  )

  # Test invalid mesh parameters
  expect_error(
    prepare_inla_data_stack(
      input_data = test_input_data,
      id_raster = test_id_rast,
      covariates = test_covariates,
      mesh_max_edge = "not numeric"
    )
  )
})

test_that("prepare_inla_data_stack handles edge cases", {
  # Test with single observation
  single_data <- test_input_data[1]
  result_single <- prepare_inla_data_stack(
    input_data = single_data,
    id_raster = test_id_rast,
    covariates = test_covariates
  )
  expect_true(is.list(result_single))

  # Test with single covariate
  single_cov <- list(cov1 = test_covariate)
  result_single_cov <- prepare_inla_data_stack(
    input_data = test_input_data,
    id_raster = test_id_rast,
    covariates = single_cov
  )
  expect_true(is.list(result_single_cov))
})
