#' Make time stamp
#'
#' @description Create a string time stamp based on current detailed date/time
#'
#' @param suffix (`character(1)`, default NULL) suffix to append to the time stamp. Useful
#'   when running batches of related models
#' @param milliseconds (`logical(1)`, default TRUE) Should milliseconds be appended to
#'   the timestamp? Useful when launching many models in quick succession.
#'
#' @return A string formatted as `'YYYYMMDD_HH_MM_SS(_optional MS)(_optional suffix)'`
#'
#' @export
make_time_stamp <- function(suffix = NULL, milliseconds = TRUE){
  if(milliseconds){
    time_stamp <- strftime(x = Sys.time(), format = '%Y%m%d_%H_%M_%OS3') |>
      gsub(pattern = '\\.', replacement = '_')
  } else {
    time_stamp <- strftime(x = Sys.time(), format = '%Y%m%d_%H_%M_%S')
  }

  # Suffix must either be NULL or length 1
  if(length(suffix) > 1) stop("suffix should be NULL or a character vector of length 1.")
  if(length(suffix) == 1) time_stamp <- paste0(time_stamp, '_', suffix)

  return(time_stamp)
}


#' Dissolve sf object by attribute
#'
#' @description Dissolve an SF object by attribute
#'
#' @details Inspired by [spatialEco::sf_dissolve]
#'
#' @param x ([sf::sf] object) SF object to dissolve
#' @param by (`character(N)`, default character(0)) Attributes to dissolve by
#'
#' @return Dissolved [sf::sf] object
#'
#' @importFrom sf st_drop_geometry
#' @export
dissolve_sf_by_attribute <- function(x, by = character(0)){
  if(length(by) == 0){
    # Dissolve all
    dissolved <- x |> sf::st_union() |> sf::st_as_sf()
  } else {
    # Dissolve by attributes
    dissolved_groups <- sf::st_drop_geometry(x)[, by] |> unique()
    dissolved <- lapply(seq_len(nrow(dissolved_groups)), function(row_id){
      dissolved_row <- dissolved_groups[row_id, ]
      sf::st_geometry(dissolved_row) <- merge(x = x, y = dissolved_row, by = by) |>
        sf::st_make_valid() |>
        sf::st_union()
      return(dissolved_row)
    }) |> do.call(what = 'rbind')
  }
  return(dissolved)
}
