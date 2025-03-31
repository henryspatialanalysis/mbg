#' Aggregate a raster to polygons: validation
#'
#' @description Data validation for aggregate_raster_to_polygons
#'
#' @seealso aggregate_raster_to_polygons
#'
#' @param data_raster [terra::SpatRaster] containing data to be aggregated to polygons.
#' @param aggregation_table Aggregation table linking pixels to polygon identifiers,
#'   created using `build_aggregation_table()`
#' @param aggregation_cols (character vector, default 'polygon_id') Polygon identifiers
#'   to use for aggregation.
#' @param method (character, default 'mean') Aggregation method: one of 'mean', 'sum',
#'   'weighted.mean', or 'weighted.sum'.
#' @param aggregated_field (character, default 'data') Name of the aggregated raster
#'   values in the output table.
#' @param z_dimension (vector, default NULL) If passing a `data_raster` with multiple
#'   layers, how should each layer be identified?
#' @param z_dimension_name (default 'z') The field name for the "z" dimension
#'   corresponding to each layer of the `data_raster`.
#' @param weighting_raster ([terra::SpatRaster], default NULL) The relative weighting of
#'   each whole pixel to the overall polygon value, for example, if calculating a
#'   population-weighted mean.
#' @param na.rm (bool, default TRUE) How to handle NA pixels in `data_raster` and
#'   `weighting_raster`.
#'
#' @return Errors if checks fail; silently passes if checks pass
#'
#' @keywords internal
#'
#' @importFrom assertthat assert_that has_name noNA
#' @importFrom terra compareGeom ncell nlyr
aggregate_raster_to_polygons_validation <- function(
  data_raster, aggregation_table, aggregation_cols, method, aggregated_field, z_dimension,
  z_dimension_name, weighting_raster, na.rm
){
  # Overload a data.table variable to pass R CMD check
  pixel_id <- NULL

  # Checks on the aggregation table
  assertthat::assert_that(inherits(aggregation_table, 'data.table'))
  assertthat::assert_that(assertthat::has_name(aggregation_table, aggregation_cols))
  # Checks on the data raster
  assertthat::assert_that(inherits(data_raster, 'SpatRaster'))
  assertthat::assert_that(assertthat::noNA(aggregation_table$pixel_id))
  max_pixel_id <- aggregation_table[, max(pixel_id)]
  assertthat::assert_that(
    (terra::ncell(data_raster)) >= max_pixel_id,
    msg = 'The data raster has dimensions larger than expected by the aggregation table.'
  )
  # Checks on the aggregation method
  valid_methods <- c('mean', 'weighted.mean', 'sum', 'weighted.sum')
  assertthat::assert_that(length(method) == 1L)
  assertthat::assert_that(method %in% valid_methods)
  # Checks on the aggregated field name
  assertthat::assert_that(inherits(aggregated_field, 'character'))
  assertthat::assert_that(length(aggregated_field) == 1L)
  # Checks on the Z dimension
  if(!is.null(z_dimension)){
    assertthat::assert_that(length(z_dimension) == terra::nlyr(data_raster))
  }
  # Check for duplicates in the aggregated table fields
  aggregated_all_fields <- c(aggregation_cols, aggregated_field)
  if(!is.null(z_dimension) | (terra::nlyr(data_raster) > 1L)){
    assertthat::assert_that(!is.null(z_dimension_name))
    aggregated_all_fields <- c(aggregated_all_fields, z_dimension_name)
  }
  duplicates <- aggregated_all_fields[duplicated(aggregated_all_fields)] |> unique()
  if(length(duplicates) > 0L) stop(
    "Fields ", paste(duplicates, collapse = ', '), " would appear multiple times in the ",
    "aggregated table."
  )
  # Checks on the weighting raster, if needed
  if(grepl('^weighted', method)){
    assertthat::assert_that(inherits(weighting_raster, 'SpatRaster'))
    assertthat::assert_that(terra::compareGeom(x = data_raster, y = weighting_raster))
    data_raster_n_layers <- terra::nlyr(data_raster)
    weighting_raster_n_layers <- terra::nlyr(weighting_raster)
    # The weighting raster will usually have the same number of layers as the weighting
    #  raster. A weighting raster with one layer will pass checks but yield a warning
    if((data_raster_n_layers > 1L) & (weighting_raster_n_layers == 1L)){
      warning(
        "Data raster has multiple layers but weighting raster has only one layer. The ",
        "same weighting raster will be used for each z dimension."
      )
    } else {
      assertthat::assert_that(data_raster_n_layers == weighting_raster_n_layers)
    }
  }
  # Checks on na.rm
  assertthat::assert_that(inherits(na.rm, 'logical'))
  assertthat::assert_that(length(na.rm) == 1L)
  invisible(NULL)
}


#' Aggregate a raster to polygons
#'
#' @description Aggregate raster values to polygons using an aggregation table
#'
#' @details This is a more efficient and feature-rich alternative to terra's zonal
#'   statistics functions. Features include:
#'   - Always fractionally aggregate, weighting by area of the pixel in a polygon
#'   - Optionally weight by both area fraction and a weighting raster (e.g. population)
#'   - Means or sums of raster values across polygons
#'   - Optionally aggregate multiple years of raster data at once
#'
#' @param data_raster [terra::SpatRaster] containing data to be aggregated to polygons.
#' @param aggregation_table [data.table::data.table] Aggregation table linking pixels to
#'   polygon identifiers, created using `build_aggregation_table()`
#' @param aggregation_cols (character vector, default 'polygon_id') Polygon identifiers
#'   to use for aggregation. Must be field names within `aggregation_table`.
#' @param method (character, default 'mean') Aggregation method: one of 'mean', 'sum',
#'   'weighted.mean', or 'weighted.sum'. The latter two methods require a
#'   `weighting_raster.`
#' @param aggregated_field (character, default 'data') Name of the aggregated raster
#'   values in the output table.
#' @param z_dimension (vector, default NULL) If passing a `data_raster` with multiple
#'   layers, how should each layer be identified? Should be a vector with length equal to
#'   the number of layers in `data_raster`. If left as `NULL`, the default, and
#'   `data_raster` has 1 layer, no z dimension will be added. If left as `NULL` and
#'   `data_raster` has more than 1 layer, will default to (1, 2, ..., N layers).
#' @param z_dimension_name (default 'z') The field name for the "z" dimension
#'   corresponding to each layer of the `data_raster`. This field is only added if
#'   `z_dimension` is passed or if `data_raster` has more than one layer.
#' @param weighting_raster ([terra::SpatRaster], default NULL) The relative weighting of
#'   each whole pixel to the overall polygon value, for example, if calculating a
#'   population-weighted mean. Required for methods 'weighted.mean' and 'weighted.sum',
#'   ignored for the other methods.
#' @param na.rm (bool, default TRUE) How to handle NA pixels in `data_raster` and
#'   `weighting_raster`. If set to TRUE but ALL pixels in a polygon are NA, will still
#'   return an NA value for the polygon.
#'
#' @return data.table containing polygon identifiers, (optionally) layer identifiers in
#'   the `z_dimension_name` column, and data values aggregated by polygon.
#'
#' @examples
#' \dontrun{
#'   polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
#'   id_raster <- build_id_raster(polygons)
#'   n_data_pixels <- sum(!is.na(terra::values(id_raster)))
#'   # Example ID raster filled with data
#'   # This is an example of pixel-level covariate data or model estimates
#'   data_raster <- mbg::values_to_raster(stats::rnorm(n_data_pixels), id_raster)
#'   # Build aggregation table, which can be used across many aggregations
#'   aggregation_table <- build_aggregation_table(
#'     polygons, id_raster, polygon_id_field = 'commune_code'
#'   )
#'   # Aggregate the raster to the polygons
#'    aggregated <- aggregate_raster_to_polygons(
#'      data_raster = data_raster,
#'      aggregation_table = aggregation_table,
#'      aggregation_cols = 'commune_code',
#'      method = 'mean'
#'    )
#'    head(aggregated)
#' }
#'
#' @seealso build_aggregation_table
#'
#' @concept aggregation
#'
#' @importFrom terra nlyr ncell values
#' @importFrom stats weighted.mean
#' @import data.table
#' @export
aggregate_raster_to_polygons <- function(
  data_raster, aggregation_table, aggregation_cols = 'polygon_id', method = 'mean',
  aggregated_field = 'data', z_dimension = NULL, z_dimension_name = 'z',
  weighting_raster = NULL, na.rm = TRUE
){
  # Overload some data.table variables to pass R CMD check
  . <- area_fraction <- w__ <- val__ <- z__ <- NULL

  ## Data preparation
  # Validate function inputs
  aggregate_raster_to_polygons_validation(
    data_raster = data_raster,
    aggregation_table = aggregation_table,
    aggregation_cols = aggregation_cols,
    method = method,
    aggregated_field = aggregated_field,
    z_dimension = z_dimension,
    z_dimension_name = z_dimension_name,
    weighting_raster = weighting_raster,
    na.rm = na.rm
  )
  # Determine whether Z dimension will be reported
  num_z <- terra::nlyr(data_raster)
  report_z <- ((num_z > 1L) | !is.null(z_dimension))
  # Fill Z dimension, if needed
  if(is.null(z_dimension)) z_dimension <- seq_len(num_z)

  ## Create a table indexing pixels to their ID
  px_table <- data.table::data.table(val__ = terra::values(data_raster, mat = F))
  px_table$pixel_id <- rep(seq_len(terra::ncell(data_raster)), times = num_z)
  px_table$z__ <- rep(z_dimension, each = terra::ncell(data_raster))
  # Optionally add weights
  need_weights <- grepl('^weighted\\.', method)
  if(need_weights){
    weights_num_z <- terra::nlyr(weighting_raster)
    if(num_z == weights_num_z){
      px_table$w__ <- terra::values(weighting_raster, mat = F)
    } else if(weights_num_z == 1){
      px_table$w__ <- rep(terra::values(weighting_raster, mat = F), times = num_z)
    } else {
      stop("ISSUE: Weighting raster has a different number of layers from data raster.")
    }
  } else {
    px_table$w__ <- 1L
  }

  ## Merge pixel table with aggregation table
  disaggregated_table <- merge(
    x = px_table,
    y = aggregation_table,
    by = 'pixel_id',
    allow.cartesian = TRUE
  )

  ## Aggregation varies depending on method
  if(method %in% c('mean', 'weighted.mean')){
    prepped_table <- disaggregated_table[
      , .(val__ = weighted.mean(val__, w = area_fraction * w__, na.rm = na.rm)),
      by = c(aggregation_cols, 'z__')
    ]
  } else if(method %in% c('sum', 'weighted.sum')){
    prepped_table <- disaggregated_table[
      , .(val__ = sum(val__ * area_fraction * w__, na.rm = na.rm)),
      by = c(aggregation_cols, 'z__')
    ]
  } else {
    stop("Method ", method, " not yet available as an aggregation function.")
  }

  ## Format for return
  if(report_z){
    if(z_dimension_name != 'z__') setnames(prepped_table, 'z__', z_dimension_name)
    setorderv(prepped_table, c(z_dimension_name, aggregation_cols))
  } else {
    prepped_table[, z__ := NULL ]
    setorderv(prepped_table, aggregation_cols)
  }
  if(aggregated_field != 'val__') setnames(prepped_table, 'val__', aggregated_field)

  # Return
  return(prepped_table)
}


#' Aggregate grid cell draws to polygons: validation
#'
#' @description Data validation for aggregate_draws_to_polygons
#'
#' @seealso aggregate_draws_to_polygons
#'
#' @param draws_matrix `matrix`, `array`, or `data.frame` corresponding to grid cell draws
#'   that will be aggregated to polygons.
#' @param aggregation_table [data.table::data.table] Aggregation table linking pixels to
#'   polygon identifiers, created using `build_aggregation_table()`
#' @param aggregation_cols (character vector, default 'polygon_id') Polygon identifiers
#'   to use for aggregation.
#' @param method (character, default 'mean') Aggregation method: one of 'mean', 'sum',
#'   'weighted.mean', or 'weighted.sum'.
#' @param z_dimension (vector, default NULL) If passing a `draws_matrix` with multiple
#'   sets of estimates, how should each layer be identified?
#' @param z_dimension_name (default 'z') The field name for the "z" dimension
#'   corresponding to each set of estimates contained in `draws_matrix`.
#' @param weighting_raster ([terra::SpatRaster], default NULL) The relative weighting of
#'   each whole pixel to the overall polygon value, for example, if calculating a
#'   population-weighted mean.
#' @param na.rm (bool, default TRUE) How to handle NA values in `draws_matrix` and
#'   `weighting_raster`.
#'
#' @return Errors if checks fail; silently passes if checks pass
#'
#' @keywords internal
#'
#' @importFrom assertthat assert_that has_name noNA
#' @importFrom terra compareGeom ncell nlyr
aggregate_draws_to_polygons_validation <- function(
  draws_matrix, aggregation_table, aggregation_cols, method, z_dimension,
  z_dimension_name, weighting_raster, na.rm
){
  # Checks on the aggregation table
  assertthat::assert_that(inherits(aggregation_table, 'data.table'))
  assertthat::assert_that(assertthat::has_name(aggregation_table, aggregation_cols))
  assertthat::assert_that(assertthat::noNA(aggregation_table$masked_pixel_id))
  non_na_px <- max(aggregation_table$masked_pixel_id)
  # Checks on the draws matrix
  assertthat::assert_that(
    inherits(draws_matrix, 'matrix') | inherits(draws_matrix, 'data.frame')
  )
  assertthat::assert_that(
    nrow(draws_matrix) %% non_na_px == 0,
    msg = paste(
      "Number of rows in the draws matrix must be exactly divisible by the number of ",
      "non-NA pixels in the ID raster, that is, max(aggregation_table$masked_pixel_id)."
    )
  )
  num_z <- nrow(draws_matrix) / non_na_px
  # Checks on the aggregation method
  valid_methods <- c('mean', 'weighted.mean', 'sum', 'weighted.sum')
  assertthat::assert_that(length(method) == 1L)
  assertthat::assert_that(method %in% valid_methods)
  # Checks on the Z dimension
  if(!is.null(z_dimension)){
    assertthat::assert_that(length(z_dimension) == num_z)
  }
  # Check for duplicates in the aggregated table fields
  num_draws <- ncol(draws_matrix)
  draw_cols <- paste0('draw_', seq_len(num_draws))
  aggregated_all_fields <- c(aggregation_cols, draw_cols)
  if(!is.null(z_dimension) | (num_z > 1L)){
    assertthat::assert_that(length(z_dimension_name) == 1L)
    aggregated_all_fields <- c(aggregated_all_fields, z_dimension_name)
  }
  duplicates <- aggregated_all_fields[duplicated(aggregated_all_fields)] |> unique()
  if(length(duplicates) > 0L) stop(
    "Fields ", paste(duplicates, collapse = ', '), " would appear multiple times in the ",
    "aggregated table."
  )
  # Checks on the weighting raster, if needed
  if(grepl('^weighted', method)){
    assertthat::assert_that(inherits(weighting_raster, 'SpatRaster'))
    weighting_raster_n_layers <- terra::nlyr(weighting_raster)
    # The number of layers in the weighting raster will usually equal the number of Z
    #  dimensions. A weighting raster with one layer will pass checks but yield a warning
    if((num_z > 1L) & (weighting_raster_n_layers == 1L)){
      warning(
        "Data raster has multiple layers but weighting raster has only one layer. The ",
        "same weighting raster will be used for each z dimension."
      )
    } else {
      assertthat::assert_that(num_z == weighting_raster_n_layers)
    }
  }
  # Checks on na.rm
  assertthat::assert_that(inherits(na.rm, 'logical'))
  assertthat::assert_that(length(na.rm) == 1L)
  invisible(NULL)
}


#' Aggregate grid cell draws to polygons
#'
#' @description Aggregate grid cell draws to polygons using an aggregation table
#'
#' @details This is a more efficient and feature-rich alternative to terra's zonal
#'   statistics functions. Features include:
#'   - Always fractionally aggregate, weighting by area of the pixel in a polygon
#'   - Optionally weight by both area fraction and a weighting raster (e.g. population)
#'   - Means or sums of raster values across polygons
#'   - Optionally aggregate multiple years of raster data at once
#'
#' @param draws_matrix `matrix`, `array`, or `data.frame` corresponding to grid cell draws
#'   that will be aggregated to polygons:\cr
#'   - Each row represents a non-NA grid cell in the ID raster. If the matrix contains
#'     multiple years of estimates, the matrix is ordered by year, then by
#'     masked_pixel_id. For example, if there are 200 non-NA pixels in the ID raster and
#'     five years of draws, then the matrix contains 1000 rows: row 200 corresponds to
#'     (year 1, masked_pixel_id 200), row 201 corresponds to (year 2, masked_pixel_id 1),
#'     and so on.\cr
#'   - Each column represents a draw. There should be no non-draw columns (such as ID
#'     fields) in the `draws_matrix`.
#' @param aggregation_table [data.table::data.table] Aggregation table linking pixels to
#'   polygon identifiers, created using `build_aggregation_table()`
#' @param aggregation_cols (character vector, default 'polygon_id') Polygon identifiers
#'   to use for aggregation. Must be field names within `aggregation_table`.
#' @param method (character, default 'mean') Aggregation method: one of 'mean', 'sum',
#'   'weighted.mean', or 'weighted.sum'. The latter two methods require a
#'   `weighting_raster.`
#' @param z_dimension (vector, default NULL) If passing a `data_raster` with multiple
#'   layers, how should each layer be identified? Should be a vector with length equal to
#'   the number of layers in `data_raster`. If left as `NULL`, the default, and
#'   `data_raster` has 1 layer, no z dimension will be added. If left as `NULL` and
#'   `data_raster` has more than 1 layer, will default to (1, 2, ..., N layers).
#' @param z_dimension_name (default 'z') The field name for the "z" dimension
#'   corresponding to each layer of the `data_raster`. This field is only added if
#'   `z_dimension` is passed or if `data_raster` has more than one layer.
#' @param weighting_raster ([terra::SpatRaster], default NULL) The relative weighting of
#'   each whole pixel to the overall polygon value, for example, if calculating a
#'   population-weighted mean. Required for methods 'weighted.mean' and 'weighted.sum',
#'   ignored for the other methods.
#' @param na.rm (bool, default TRUE) How to handle NA pixels in `data_raster` and
#'   `weighting_raster`. If set to TRUE but ALL pixels in a polygon are NA, will still
#'   return an NA value for the polygon.
#'
#' @return [data.table::data.table] containing polygon identifiers, (optionally) layer
#'   identifiers in the `z_dimension_name` column, and data values aggregated by polygon.
#'
#' @examples
#' \dontrun{
#'   polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
#'   id_raster <- build_id_raster(polygons)
#'   n_data_pixels <- sum(!is.na(terra::values(id_raster)))
#'   # Example grid-level draws from e.g. mbg::generate_cell_draws_and_summarize()
#'   draws_matrix <- matrix(rnorm(n_data_pixels * 5), nrow = n_data_pixels)
#'   # Build aggregation table, which can be used across many aggregations
#'   aggregation_table <- build_aggregation_table(
#'     polygons, id_raster, polygon_id_field = 'commune_code'
#'   )
#'   # Aggregate the grid-level draws to polygon-level draws
#'   aggregated <- aggregate_draws_to_polygons(
#'     draws_matrix = draws_matrix,
#'     aggregation_table = aggregation_table,
#'     aggregation_cols = 'commune_code',
#'     method = 'mean'
#'   )
#'   head(aggregated)
#' }
#'
#' @concept aggregation
#'
#' @seealso build_aggregation_table
#'
#' @importFrom terra nlyr ncell values
#' @importFrom stats na.omit weighted.mean
#' @import data.table
#' @export
aggregate_draws_to_polygons <- function(
  draws_matrix, aggregation_table, aggregation_cols = 'polygon_id', method = 'mean',
  z_dimension = NULL, z_dimension_name = 'z', weighting_raster = NULL, na.rm = TRUE
){
  # Overload some data.table variables to pass R CMD check
  . <- pixel_id <- masked_pixel_id <- i.masked_pixel_id <- w__ <- i.w__ <-
    area_fraction <- z__ <- NULL

  ## Data preparation
  # Validate function inputs
  aggregate_draws_to_polygons_validation(
    draws_matrix = draws_matrix,
    aggregation_table = aggregation_table,
    aggregation_cols = aggregation_cols,
    method = method,
    z_dimension = z_dimension,
    z_dimension_name = z_dimension_name,
    weighting_raster = weighting_raster,
    na.rm = na.rm
  )
  # Determine whether Z dimension will be reported
  num_non_na_pixels <- max(aggregation_table$masked_pixel_id)
  num_z <- nrow(draws_matrix) / num_non_na_pixels
  report_z <- ((num_z > 1L) | !is.null(z_dimension))
  # Fill Z dimension, if needed
  if(is.null(z_dimension)) z_dimension <- seq_len(num_z)

  ## Create a table indexing draws to a masked pixel ID
  draws_table <- as.data.table(draws_matrix)
  num_draws <- ncol(draws_table)
  draw_cols <- paste0('draw_', seq_len(num_draws))
  colnames(draws_table) <- draw_cols
  # Add masked pixel ID
  draws_table$masked_pixel_id <- rep(seq_len(num_non_na_pixels), times = num_z)
  draws_table$z__ <- rep(z_dimension, each = num_non_na_pixels)
  # Optionally add weights
  need_weights <- grepl('^weighted\\.', method)
  if(need_weights){
    weights_num_z <- terra::nlyr(weighting_raster)
    weights_ncell <- terra::ncell(weighting_raster)
    id_to_masked_id <- unique(aggregation_table[, .(pixel_id, masked_pixel_id)])
    weights_table <- (
      data.table(
        w__ = terra::values(weighting_raster, mat = F),
        pixel_id = rep(seq_len(weights_ncell), times = weights_num_z),
        z__ = rep(z_dimension[seq_len(weights_num_z)], each = weights_ncell)
      )
      [id_to_masked_id, masked_pixel_id := i.masked_pixel_id, on = 'pixel_id']
      [!is.na(masked_pixel_id), ]
    )
    if(weights_num_z == num_z){
      draws_table[weights_table, w__ := i.w__, on = c('masked_pixel_id', 'z__')]
    } else if(weights_num_z == 1){
      draws_table[weights_table, w__ := i.w__, on = 'masked_pixel_id']
    } else {
      stop("ISSUE: Weighting raster has a different number of layers from data raster.")
    }
  } else {
    draws_table$w__ <- 1L
  }

  ## Merge pixel table with aggregation table
  disaggregated_table <- merge(
    x = draws_table,
    y = aggregation_table,
    by = 'masked_pixel_id',
    allow.cartesian = T
  )

  ## Aggregation varies depending on method
  if(method %in% c('mean', 'weighted.mean')){
    prepped_table <- disaggregated_table[
      , lapply(.SD, weighted.mean, w = area_fraction * w__, na.rm = na.rm),
      .SDcols = draw_cols,
      by = c(aggregation_cols, 'z__')
    ]
  } else if(method %in% c('sum', 'weighted.sum')){
    prepped_table <- disaggregated_table[
      , lapply(.SD, function(x) sum(x * area_fraction * w__, na.rm = na.rm)),
      .SDcols = draw_cols,
      by = c(aggregation_cols, 'z__')
    ]
  } else {
    stop("Method ", method, " not yet available as an aggregation function.")
  }

  ## Format for return
  if(report_z){
    if(z_dimension_name != 'z__') setnames(prepped_table, 'z__', z_dimension_name)
    setorderv(prepped_table, c(z_dimension_name, aggregation_cols))
  } else {
    prepped_table[, z__ := NULL ]
    setorderv(prepped_table, aggregation_cols)
  }

  # Return
  return(prepped_table)
}
