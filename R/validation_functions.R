#' Generate RMSE from an estimated raster surface and point data
#'
#' @details For examples, see \code{vignette('model-comparison', package = 'mbg')}
#'
#' @param estimates ([terra::SpatRaster]) Raster surface containing point estimates. This
#'   could also be the mean surface of a Bayesian geostatistical model
#' @param validation_data (`data.frame`)\cr
#'   Table containing at least the following fields:\cr
#'   * x (`numeric`) location x position, in the same projection as `estimates`\cr
#'   * y (`numeric`) location y position, in the same projection as `estimates`\cr
#'   * (Outcome field) See below
#' @param outcome_field (`character(1)`) Column in `validation_data` containing the values
#'   that should be compared against the `estimates` raster surface.
#' @param na.rm (`logical(1)`, default FALSE) Should NA values be dropped from the RMSE
#'   calculation?
#'
#' @concept validation
#'
#' @return A single number giving RMSE between the point data and estimates raster.
#' @importFrom terra extract
#' @export
rmse_raster_to_point <- function(
  estimates,
  validation_data,
  outcome_field,
  na.rm = FALSE
){
  point_estimates <- terra::extract(
    x = estimates,
    y = validation_data[, c('x', 'y')],
    ID = FALSE
  ) |> unlist()
  rmse <- (point_estimates - validation_data[[outcome_field]])**2 |>
    mean(na.rm = na.rm) |>
    sqrt()
  return(rmse)
}


#' Generate log posterior predictive density from a geostatistical surface onto point data
#'
#' @details Calculated across draws. Requires an ID raster to match each point observation
#'   to a set of draws. Assumes binomial data.
#'
#' @details For examples, see \code{vignette('model-comparison', package = 'mbg')}
#'
#' @param draws (`matrix`) A predictive draw matrix, where each row corresponds to a
#'   pixel in the `id_raster` and each column corresponds to one sampled estimate of the
#'   outcome.
#' @param validation_data (`data.frame`) Table containing at least the following
#'  fields:\cr
#'   * x (`numeric`) location x position, in the same projection as `id_raster`\cr
#'   * y (`numeric`) location y position, in the same projection as `id_raster`\cr
#'   * indicator (`integer`) The number of events in the population\cr
#'   * samplesize (`integer`) The total population, denominator for `indicator`
#' @param id_raster ([terra::SpatRaster]) Raster showing the sample study area, created
#'   using [build_id_raster].
#' @param na.rm (`logical(1)`, default FALSE) Should NA values be omitted from the LPD
#'   calculation?
#'
#' @return (`numeric(1)`) Log predictive density of the validation data given the draw
#'   estimates.
#'
#' @concept validation
#'
#' @import data.table
#' @importFrom terra extract
#' @importFrom stats dbinom
#' @export
log_posterior_density <- function(draws, validation_data, id_raster, na.rm = FALSE){
  # Overload some data.table variables to pass R CMD check
  draws_row <- i.draws_row <- NULL

  # Get the row of predictive draws corresponding to each data point
  pixel_ids <- terra::extract(
    x = id_raster,
    y = validation_data[, c('x', 'y')],
    ID = FALSE,
    dataframe = TRUE
  ) |>
    data.table::as.data.table() |>
    data.table::setnames(new = 'pixel_id')
  row_lookup <- terra::values(
    x = id_raster,
    dataframe = TRUE,
    na.rm = TRUE
  ) |>
    data.table::as.data.table() |>
    data.table::setnames(new = 'pixel_id')
  row_lookup[, draws_row := .I ]
  pixel_ids[row_lookup, draws_row := i.draws_row, on = 'pixel_id']
  # Get log density at each data point
  log_density_by_observation <- sapply(seq_len(nrow(validation_data)), function(row_ii){
    stats::dbinom(
      x = validation_data$indicator[row_ii],
      size = validation_data$samplesize[row_ii],
      prob = draws[pixel_ids$draws_row[row_ii], ]
    ) |>
      mean(na.rm = na.rm) |>
      log()
  })
  log_density <- sum(log_density_by_observation, na.rm = na.rm)
  return(log_density)
}
