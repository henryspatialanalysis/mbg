# Aggregate grid cell draws to polygons

Aggregate grid cell draws to polygons using an aggregation table

## Usage

``` r
aggregate_draws_to_polygons(
  draws_matrix,
  aggregation_table,
  aggregation_cols = "polygon_id",
  method = "mean",
  z_dimension = NULL,
  z_dimension_name = "z",
  weighting_raster = NULL,
  na.rm = TRUE
)
```

## Arguments

- draws_matrix:

  `matrix`, `array`, or `data.frame` corresponding to grid cell draws
  that will be aggregated to polygons:  

  - Each row represents a non-NA grid cell in the ID raster. If the
    matrix contains multiple years of estimates, the matrix is ordered
    by year, then by masked_pixel_id. For example, if there are 200
    non-NA pixels in the ID raster and five years of draws, then the
    matrix contains 1000 rows: row 200 corresponds to (year 1,
    masked_pixel_id 200), row 201 corresponds to (year 2,
    masked_pixel_id 1), and so on.  

  - Each column represents a draw. There should be no non-draw columns
    (such as ID fields) in the `draws_matrix`.

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

[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
containing polygon identifiers, (optionally) layer identifiers in the
`z_dimension_name` column, and data values aggregated by polygon.

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
  # Example grid-level draws from e.g. mbg::generate_cell_draws_and_summarize()
  draws_matrix <- matrix(rnorm(n_data_pixels * 5), nrow = n_data_pixels)
  # Build aggregation table, which can be used across many aggregations
  aggregation_table <- build_aggregation_table(
    polygons, id_raster, polygon_id_field = 'commune_code'
  )
  # Aggregate the grid-level draws to polygon-level draws
  aggregated <- aggregate_draws_to_polygons(
    draws_matrix = draws_matrix,
    aggregation_table = aggregation_table,
    aggregation_cols = 'commune_code',
    method = 'mean'
  )
  head(aggregated)
} # }
```
