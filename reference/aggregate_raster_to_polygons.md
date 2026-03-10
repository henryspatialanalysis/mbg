# Aggregate a raster to polygons

Aggregate raster values to polygons using an aggregation table

## Usage

``` r
aggregate_raster_to_polygons(
  data_raster,
  aggregation_table,
  aggregation_cols = "polygon_id",
  method = "mean",
  aggregated_field = "data",
  z_dimension = NULL,
  z_dimension_name = "z",
  weighting_raster = NULL,
  na.rm = TRUE
)
```

## Arguments

- data_raster:

  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  containing data to be aggregated to polygons.

- aggregation_table:

  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  Aggregation table linking pixels to polygon identifiers, created using
  [`build_aggregation_table()`](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md)

- aggregation_cols:

  (character vector, default 'polygon_id') Polygon identifiers to use
  for aggregation. Must be field names within `aggregation_table`.

- method:

  (character, default 'mean') Aggregation method: one of 'mean', 'sum',
  'weighted.mean', or 'weighted.sum'. The latter two methods require a
  `weighting_raster.`

- aggregated_field:

  (character, default 'data') Name of the aggregated raster values in
  the output table.

- z_dimension:

  (vector, default NULL) If passing a `data_raster` with multiple
  layers, how should each layer be identified? Should be a vector with
  length equal to the number of layers in `data_raster`. If left as
  `NULL`, the default, and `data_raster` has 1 layer, no z dimension
  will be added. If left as `NULL` and `data_raster` has more than 1
  layer, will default to (1, 2, ..., N layers).

- z_dimension_name:

  (default 'z') The field name for the "z" dimension corresponding to
  each layer of the `data_raster`. This field is only added if
  `z_dimension` is passed or if `data_raster` has more than one layer.

- weighting_raster:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  default NULL) The relative weighting of each whole pixel to the
  overall polygon value, for example, if calculating a
  population-weighted mean. Required for methods 'weighted.mean' and
  'weighted.sum', ignored for the other methods.

- na.rm:

  (bool, default TRUE) How to handle NA pixels in `data_raster` and
  `weighting_raster`. If set to TRUE but ALL pixels in a polygon are NA,
  will still return an NA value for the polygon.

## Value

data.table containing polygon identifiers, (optionally) layer
identifiers in the `z_dimension_name` column, and data values aggregated
by polygon.

## Details

This is a more efficient and feature-rich alternative to terra's zonal
statistics functions. Features include:

- Always fractionally aggregate, weighting by area of the pixel in a
  polygon

- Optionally weight by both area fraction and a weighting raster (e.g.
  population)

- Means or sums of raster values across polygons

- Optionally aggregate multiple years of raster data at once

## See also

build_aggregation_table

## Examples

``` r
if (FALSE) { # \dontrun{
  polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
  id_raster <- build_id_raster(polygons)
  n_data_pixels <- sum(!is.na(terra::values(id_raster)))
  # Example ID raster filled with data
  # This is an example of pixel-level covariate data or model estimates
  data_raster <- mbg::values_to_raster(stats::rnorm(n_data_pixels), id_raster)
  # Build aggregation table, which can be used across many aggregations
  aggregation_table <- build_aggregation_table(
    polygons, id_raster, polygon_id_field = 'commune_code'
  )
  # Aggregate the raster to the polygons
   aggregated <- aggregate_raster_to_polygons(
     data_raster = data_raster,
     aggregation_table = aggregation_table,
     aggregation_cols = 'commune_code',
     method = 'mean'
   )
   head(aggregated)
} # }
```
