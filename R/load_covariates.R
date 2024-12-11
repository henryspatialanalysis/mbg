#' Helper function to get the current year
#'
#' @keywords internal
#'
get_current_year <- function(){
  return(format(Sys.Date(), '%Y') |> as.integer())
}

#' Load covariates
#'
#' @details Load a list of covariates from a specified directory structure
#'
#' @param directory Directory containing all covariate sub-directories
#' @param covariates_table `data.frame` containing at least the following fields:
#'   - 'covariate': (character): Name of the covariate
#'   - 'annual': (logical) Does the covariate vary by year? If so, look for the `year` in
#'     the name of the file.
#'   - 'transform': (character) Name of a function to use to transform the covariate.
#'     Common options include 'identity' (no transformation), 'sqrt', 'abs', and 'log1p'
#'   - 'normalize': (logical) Should the covariate be rescaled to have mean 0 and standard
#'     deviation 1 across all pixels in the study area? Generally should be set to TRUE
#'     for predictive covariates.
#' @param id_raster [terra::SpatRaster] with non-NA pixels delineating the extent of the
#'   study area
#' @param year (`numeric`, default NULL) Year of data to for time-varying covariates.
#'   If NULL, the default, uses the current year.
#' @param file_format (`character`, default 'tif') File format for the raster covariate
#'   data. Used to search for the input file within the proper containing folder.
#' @param add_intercept (`logical`, default FALSE) Should a covariate called "intercept",
#'   a raster object with 1s in all required cells, be placed at the start of the returned
#'   covariates list?
#' @param check_previous_years (`integer` > 0, default 10) If annual data is not found in
#'   this year, how many previous years should be checked? If 0, will not check any
#'   previous years.
#'
#' @return A named list of formatted covariates. Each list item is a [terra::SpatRaster]
#'   with one layer and the same dimensions as the `id_raster`
#'
#' @concept core_inputs
#'
#' @importFrom assertthat assert_that
#' @importFrom stats sd
#' @importFrom terra rast crop mask global
#' @export
load_covariates <- function(
  directory, covariates_table, id_raster, year = NULL, file_format = 'tif',
  add_intercept = FALSE, check_previous_years = 10
){
  # Input data checks
  assertthat::assert_that(dir.exists(directory))
  assertthat::assert_that(inherits(covariates_table, 'data.frame'))
  assertthat::assert_that(
    all(c('covariate', 'annual', 'transform', 'normalize') %in% colnames(covariates_table))
  )
  assertthat::assert_that(inherits(id_raster, 'SpatRaster'))
  if(is.null(year)) year <- get_current_year()
  year <- as.integer(year)
  assertthat::assert_that(length(year) == 1)
  assertthat::assert_that(!is.na(year))
  assertthat::assert_that(length(file_format) == 1)
  assertthat::assert_that(inherits(file_format, 'character'))
  assertthat::assert_that(inherits(add_intercept, 'logical'))
  assertthat::assert_that(length(add_intercept) == 1)

  # Load each covariate based on its settings
  cov_names <- covariates_table$covariate
  covariates_list <- vector('list', length = length(cov_names))
  names(covariates_list) <- cov_names

  for(cov_name in cov_names){
    cov_settings <- as.list(covariates_table[covariates_table$covariate == cov_name, ])
    # Check settings for this covariate
    assertthat::assert_that(inherits(cov_settings$annual, 'logical'))
    assertthat::assert_that(inherits(cov_settings$transform, 'character'))
    assertthat::assert_that(inherits(cov_settings$normalize, 'logical'))
    # The covariate sub-directory is based on the covariate name
    cov_sub_dir <- file.path(directory, cov_name)
    assertthat::assert_that(dir.exists(cov_sub_dir))
    # Search for all files matching the extension in the covariate sub-directory
    file_format <- gsub("[[:punct:]]", "", file_format)
    paths <- list.files(
      cov_sub_dir,
      pattern = paste0("\\.", file_format, "$"),
      ignore.case = TRUE
    )
    # For annual covariates, the year must also be in the filename
    if(cov_settings$annual){
      test_year <- year
      annual_paths <- grep(pattern = paste0('_', test_year, '_'), x = paths, value = TRUE)
      # Check this year (or previous years, if specified) for annual data
      while((length(annual_paths) == 0) & (test_year > (year - check_previous_years))){
        test_year <- test_year - 1
        annual_paths <- grep(pattern = paste0('_', test_year, '_'), x = paths, value = TRUE)
      }
      if((length(annual_paths) > 0) & (test_year != year)) warning(paste(
        "Could not find", year, "data for covariate", cov_name, "-- using", test_year,
        "data instead."
      ))
      paths <- annual_paths
    }
    # The search must match exactly one file
    if(length(paths) == 0) stop("File search for covariate ", cov_name, " matched no files.")
    if(length(paths) > 1) stop(
      "File search for covariate ", cov_name, " matched multiple files:\n",
      paste(paths, collapse='\n')
    )
    # Otherwise, load the covariate
    this_cov <- (
      file.path(cov_sub_dir, paths[1]) |>
      terra::rast(win = terra::ext(id_raster), snap = 'near') |>
      terra::crop(y = id_raster) |>
      terra::mask(mask = id_raster)
    )
    assertthat::assert_that(terra::nlyr(this_cov) == 1L)
    # Optionally transform the covariate
    this_cov <- get(cov_settings$transform)(this_cov)
    # Optionally normalize the covariate
    if(cov_settings$normalize){
      this_cov <- (
        (this_cov - terra::global(this_cov, mean, na.rm = T)[1, 1]) /
        terra::global(this_cov, sd, na.rm = T)[1, 1]
      )
    }
    # If the covariate loads successfully and passes all checks, add it to the list
    covariates_list[[cov_name]] <- this_cov
  }

  # Optionally add the intercept
  if(add_intercept){
    covariates_list <- c(list(intercept = id_raster * 0 + 1), covariates_list)
  }

  # Return
  return(covariates_list)
}
