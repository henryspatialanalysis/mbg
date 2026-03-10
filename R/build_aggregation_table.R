
#' Validation: Build aggregation table
#'
#' @description Input data validation for `build_aggregation_table()`
#'
#' @param polygons [terra::SpatVector] object. Should contain a unique ID field.
#' @param id_raster [terra::SpatRaster] object. ID raster created by `build_id_raster()`
#'   for the polygons object. Should have the same CRS as `polygons` and completely cover
#'   it.
#' @param polygon_id_field (`character`) Unique identifier field in `polygons`.
#' @param polygon_ids (vector, any type) Polygon identifiers from `polygon_id_field`.
#'
#' @return Errors if checks fail; silently passes if checks pass
#'
#' @seealso build_aggregation_table
#'
#' @keywords internal
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
#' @param polygons [terra::SpatVector] object. Should contain a unique ID field.
#' @param id_raster [terra::SpatRaster] object. ID raster created by `build_id_raster()`
#'   for the polygons object. Should have the same CRS as `polygons` and completely cover
#'   it.
#' @param polygon_id_field (`character(1)`) Unique identifier field in `polygons`.
#' @param verbose (`logical(1)`, default FALSE) Show progress for building aggregation
#'   rows for each polygon?
#'
#' @return data.table with fields:
#'   * polygon_id: Unique polygon identifier
#'   * pixel_id: unique pixel ID from the ID raster
#'   * masked_pixel_id: Index counting only non-NA pixels from the ID raster
#'   * area_fraction: fraction of the pixel area falling within this polygon
#'   * Merged fields from the table of polygons
#'
#' @examples
#' \dontrun{
#'   polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
#'   id_raster <- build_id_raster(polygons)
#'   aggregation_table <- build_aggregation_table(
#'     polygons, id_raster, polygon_id_field = 'commune_code'
#'   )
#' }
#'
#' @seealso calculate_pixel_fractions_single_polygon()
#'
#' @concept aggregation
#'
#' @import data.table
#' @importFrom terra vect crop values
#' @export
build_aggregation_table <- function(
  polygons, id_raster, polygon_id_field, verbose = FALSE
){
  # Overload some data.table variables to pass R CMD check
  . <- dummy_id <- pixel_id <- masked_pixel_id <- i.masked_pixel_id <- polygon_id  <-
    area_fraction <- NULL

  # If polygons inherit sf, convert to SpatVector
  if(inherits(polygons, 'sf')) polygons <- terra::vect(polygons)
  # Input data validation
  polys_dt <- data.table::as.data.table(polygons)
  poly_ids <- polys_dt[[polygon_id_field]]
  build_aggregation_table_validation(
    polygons = polygons,
    id_raster = id_raster,
    polygon_id_field = polygon_id_field,
    polygon_ids = poly_ids
  )

  # Crop the polygons to the id raster
  polygons_cropped <- terra::crop(x = polygons, y = id_raster, ext = TRUE)
  if(verbose){
    dropped_rows <- nrow(polygons) - nrow(polygons_cropped)
    if(dropped_rows > 0L) message(paste(
      "Dropped", dropped_rows, "polygons not in the id_raster extent."
    ))
  }
  poly_ids <- as.data.frame(polygons_cropped)[[polygon_id_field]]

  agg_table <- terra::extract(
    x = id_raster,
    y = polygons_cropped,
    fun = 'table',
    exact = TRUE,
    ID = TRUE
  ) |>
    data.table::as.data.table() |>
    data.table::setnames(c('dummy_id', 'pixel_id', 'area_fraction'))
  (agg_table
    [, polygon_id := sapply(dummy_id, function(x) poly_ids[x])]
    [, dummy_id := NULL]
  )

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
  polys_dt_sub <- polys_dt[, keep_cols, with = FALSE]
  agg_table_full <- merge(
    x = agg_table[, .(polygon_id, pixel_id, masked_pixel_id, area_fraction)],
    y = polys_dt_sub,
    by = 'polygon_id'
  )[order(polygon_id, pixel_id)]

  # Return the aggregation table
  return(agg_table_full)
}
