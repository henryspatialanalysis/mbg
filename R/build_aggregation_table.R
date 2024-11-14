#' Calculate fractional pixels in a polygon
#'
#' @description Calculate the fraction of each pixel's area that falls within a single
#'   polygon
#'
#' @details This is a helper function called by `build_aggregation_table()`.
#'
#' @param polygon terra SpatVector object of length 1. The polygon to calculate fractional
#'   areas across.
#' @param id_raster terra SpatRaster object. ID raster created for the set of all polygons
#'   to be considered, created by `build_id_raster()`.
#' @param polygon_id (optional). ID for this polygon. Must have length 1.
#'
#' @return data.table containing two or three columns:
#'   * pixel_id: unique pixel ID from the ID raster
#'   * area_fraction: fraction of the pixel area falling within this polygon
#'   * polygon_id (optional): If `polygon_id` was defined, it is added to the table
#'
#' @importFrom assertthat assert_that is.scalar
#' @importFrom data.table data.table
#' @importFrom terra same.crs extract
#' @export
calculate_pixel_fractions_single_polygon <- function(polygon, id_raster, polygon_id = NULL){
  # Input data checks
  assertthat::assert_that(inherits(polygon, 'SpatVector'))
  assertthat::assert_that(inherits(id_raster, 'SpatRaster'))
  assertthat::assert_that(terra::same.crs(id_raster, polygon))
  assertthat::assert_that(nrow(polygon) == 1)
  if(!is.null(polygon_id)){
    assertthat::assert_that(assertthat::is.scalar(polygon_id))
  }

  # Use terra::extract to get pixel area fractions falling within the polygon
  # Drop the first column, which gives the polygon row ID
  extract_wide <- terra::extract(
    x = id_raster,
    y = polygon,
    fun = 'table',
    exact = TRUE,
    ID = FALSE
  )
  # Reshape long
  extract_long <- data.table::data.table(
    pixel_id = colnames(extract_wide) |> as.integer(),
    area_fraction = unlist(extract_wide[1, ])
  )
  # Optionally add the polygon ID to the table
  if(!is.null(polygon_id)) extract_long$polygon_id <- polygon_id

  return(extract_long)
}


#' Validation: Build aggregation table
#'
#' @description Input data validation for `build_aggregation_table()`
#'
#' @param polygons terra SpatVector object. Should contain a unique ID field.
#' @param id_raster terra SpatRaster object. ID raster created by `build_id_raster()` for
#'   the polygons object. Should have the same CRS as `polygons` and completely cover it.
#' @param polygon_id_field (character) Unique identifier field in `polygons`.
#' @param polygon_ids (vector, any type) Polygon identifiers from `polygon_id_field`.
#'
#' @return Errors if checks fail; silently passes if checks pass
#'
#' @seealso build_aggregation_table
#'
#' @importFrom assertthat assert_that has_name noNA
#' @importFrom terra same.crs
#' @importFrom glue glue
build_aggregation_table_validation <- function(
  polygons, id_raster, polygon_id_field, polygon_ids
){
  # Checks on the polygons spatVector
  assertthat::assert_that(inherits(polygons, 'SpatVector'))
  assertthat::assert_that(nrow(polygons) >= 1)
  # Checks on the id_raster spatRaster
  assertthat::assert_that(inherits(id_raster, 'SpatRaster'))
  assertthat::assert_that(terra::same.crs(id_raster, polygons))
  # Checks on the polygon ID field
  assertthat::assert_that(assertthat::has_name(polygons, polygon_id_field))
  # Checks on the vector of polygon IDs
  assertthat::assert_that(length(polygon_ids) == nrow(polygons))
  assertthat::assert_that(assertthat::noNA(polygon_ids))
  assertthat::assert_that(
    sum(duplicated(polygon_ids)) == 0L,
    msg = glue::glue("Polygon ID field {polygon_id_field} contains duplicates.")
  )
  invisible(NULL)
}


#' Build aggregation table
#'
#' @description Build a table to quickly aggregate from pixels to polygons
#'
#' @param polygons terra SpatVector object. Should contain a unique ID field.
#' @param id_raster terra SpatRaster object. ID raster created by `build_id_raster()` for
#'   the polygons object. Should have the same CRS as `polygons` and completely cover it.
#' @param polygon_id_field (character) Unique identifier field in `polygons`.
#' @param verbose (logical) Show progress for building aggregation rows for each polygon?
#'
#' @return data.table with fields:
#'   * polygon_id: Unique polygon identifier
#'   * pixel_id: unique pixel ID from the ID raster
#'   * masked_pixel_id: Index counting only non-NA pixels from the ID raster
#'   * area_fraction: fraction of the pixel area falling within this polygon
#'   * Merged fields from the table of polygons
#'
#' @seealso calculate_pixel_fractions_single_polygon()
#'
#' @import data.table
#' @export
build_aggregation_table <- function(polygons, id_raster, polygon_id_field, verbose = TRUE){

  polys_dt <- as.data.table(polygons)
  poly_ids <- polys_dt[[polygon_id_field]]

  # Input data validation
  build_aggregation_table_validation(
    polygons = polygons,
    id_raster = id_raster,
    polygon_id_field = polygon_id_field,
    polygon_ids = poly_ids
  )

  # Crop the polygons to the id raster
  polygons_cropped <- terra::crop(x = polygons, y = id_raster, ext = TRUE)
  poly_ids <- polys_dt[[polygon_id_field]]

  # Build the aggregation table by calling zonal statistics for each polygon
  agg_table <- lapply(poly_ids, function(poly_id){
    # Create smaller spatial objects for calculating fractional zonal statistics
    if(verbose) message('.', appendLF = F)
    one_poly <- polygons_cropped[poly_ids == poly_id, ]
    id_raster_sub <- terra::crop(x = id_raster, y = one_poly, ext = TRUE, snap = 'out')
    # Mask missing values with -1
    terra::values(id_raster_sub)[which(is.na(terra::values(id_raster_sub, mat = F)))] <- -1
    # Get fractional pixel areas for the polygon
    pixel_fractions <- calculate_pixel_fractions_single_polygon(
      polygon = one_poly,
      id_raster = id_raster_sub,
      polygon_id = poly_id
    )
    # Drop -1 (masked) pixel IDs and return
    return(pixel_fractions[pixel_id >= 0, ])
  }) |> data.table::rbindlist()
  # Add a newline to finish off the progress bar, if needed
  if(verbose) message("")

  # Merge on 'masked_pixel_id'
  masked_pixel_table <- data.table::data.table(
    pixel_id = terra::values(id_raster, mat = F, dataframe = F)
  ) |> na.omit()
  masked_pixel_table[, masked_pixel_id := .I ]
  agg_table[masked_pixel_table, masked_pixel_id := i.masked_pixel_id, on = 'pixel_id']

  # Merge on the rest of the polygon identifiers
  if(polygon_id_field != 'polygon_id') polys_dt[, polygon_id := get(polygon_id_field)]
  keep_cols <- setdiff(
    x = colnames(polys_dt),
    y = c('pixel_id', 'masked_pixel_id', 'area_fraction')
  )
  polys_dt_sub <- polys_dt[, ..keep_cols]
  agg_table_full <- merge(
    x = agg_table[, .(polygon_id, pixel_id, masked_pixel_id, area_fraction)],
    y = polys_dt_sub,
    by = 'polygon_id'
  )[order(polygon_id, pixel_id)]

  # Return the aggregation table
  return(agg_table_full)
}
