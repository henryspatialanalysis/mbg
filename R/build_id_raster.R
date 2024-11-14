#' Make world template raster
#'
#' @description Create a template raster for the world with approximately 5x5km resolution
#'   at the equator, matching many common raster covariates for health.
#'
#' @details The raster has the following specifications:
#'   * 4320 rows, 8640 columns
#'   * Resolution: 0.04166667 decimal degrees, approx. 5km at the equator
#'   * CRS: WGS 1984 unprojected latitude/longitude, `terra::crs('EPSG:4326')`
#'   * Values: All NA. Used exclusively for creating a shapefile-specific ID raster
#'
#' @returns terra SpatRaster matching the details above
#'
#' @importFrom terra rast crs
#' @export
make_world_template_raster <- function(){
  world_template_raster <- terra::rast(
    nrows = 4320,
    ncols = 8640,
    crs = terra::crs("EPSG:4326")
  )
  return(world_template_raster)
}

#' Build ID raster
#'
#' @description Build an ID raster matching the extent of a vector dataset
#'
#' @details The ID raster will be used to build the aggregation table. Each pixel has a
#'   unique integer value from 1 to the number of pixels in the ID raster.
#'
#' @param polygons terra SpatVector object. The polygons to be aggregated to
#' @param template_raster (optional) terra SpatRaster object. The template raster should
#'   contain and have the same CRS as the polygons. If template raster is `NULL`, the
#'   default, uses the default world template raster from `make_world_template_raster()`.
#'
#' @return ID raster. A terra SpatRaster object that minimally encloses the polygons
#'
#' @importFrom terra same.crs crop values
#' @importFrom assertthat assert_that
#' @export
build_id_raster <- function(polygons, template_raster = NULL){
  # Get default template raster, if needed
  if(is.null(template_raster)) template_raster <- make_world_template_raster()

  # Input data checks
  assertthat::assert_that(inherits(template_raster, 'SpatRaster'))
  assertthat::assert_that(inherits(polygons, 'SpatVector'))
  assertthat::assert_that(terra::same.crs(template_raster, polygons))

  # Set ID raster extent
  id_raster <- terra::crop(x = template_raster, y = polygons, snap = 'out')
  # Fill values of the ID raster
  terra::values(id_raster) <- seq_len(prod(dim(id_raster)))
  # Drop pixels that do not overlap the polygons
  id_raster <- terra::mask(
    x = id_raster, mask = polygons, touches = T, updatevalue = NA_integer_
  )

  return(id_raster)
}

#' Insert values into a raster
#'
#' @description Insert a vector or matrix of values into an ID spatRaster
#'
#' @details The length of the vector or matrix must be a multiple of the number of non-NA
#'   pixels in the ID raster. Values from the vector/matrix are then inserted into the
#'   non-NA pixels of the spatRaster.
#'
#' @param x Vector, matrix, data.frame, or data.table of values that will be inserted into
#'   the ID raster. The length of x must be exactly divisible by
#'   `sum(!is.na(terra::values(id_raster)))`. Data.frames are converted to matrices, and
#'   then matrices are converted to vectors using [as.matrix()] and [as.vector()]
#'   respectively before processing. For that reason, data.frames should only contain
#'   fields with values to be inserted (such as a data.frame of draws).
#' @param id_raster ID raster showing the outline of the study area, created using
#'  [build_id_raster()]. Should have 1 layer.
#'
#' @return SpatRaster with the same outline as the ID raster and (# values / # non-NA
#'   pixels in the ID raster) layers.
#'
#' @seealso [build_id_raster()]
#'
#' @importFrom terra values
#' @importFrom assertthat assert_that
#' @export
values_to_raster <- function(x, id_raster){
  # Convert x to vector
  if(inherits(x, 'data.frame')) x <- as.matrix(x)
  if(inherits(x, 'matrix')) x <- as.vector(x)

  # Input data checks
  assertthat::assert_that(
    is.vector(x),
    msg = 'x must be a vector, matrix, data.frame, or data.table of values.'
  )
  assertthat::assert_that(inherits(id_raster, 'SpatRaster'))
  assertthat::assert_that(terra::nlyr(id_raster) == 1)
  which_not_na <- which(!is.na(terra::values(id_raster)))
  n_not_na <- length(which_not_na)
  n_vals <- length(x)
  assertthat::assert_that(n_vals > 0)
  assertthat::assert_that(
    n_vals %% n_not_na == 0,
    msg = paste(
      "The length of x must be exactly divisible by the number of non-NA pixels in the",
      "ID raster."
    )
  )

  # Insert values into the non-NA pixels of the ID raster
  n_layers <- n_vals / n_not_na
  filled_rast <- lapply(seq_len(n_layers), function(layer_i){
    fill_vals <- x[seq((layer_i - 1) * n_not_na + 1, layer_i * n_not_na)]
    new_raster <- id_raster
    terra::values(new_raster)[which_not_na] <- fill_vals
    return(new_raster)
  }) |> terra::rast()

  return(filled_rast)
}
