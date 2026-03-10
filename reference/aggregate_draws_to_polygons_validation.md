# Aggregate grid cell draws to polygons: validation

Data validation for aggregate_draws_to_polygons

## Usage

``` r
aggregate_draws_to_polygons_validation(
  draws_matrix,
  aggregation_table,
  aggregation_cols,
  method,
  z_dimension,
  z_dimension_name,
  weighting_raster,
  na.rm
)
```

## Arguments

- draws_matrix:

  `matrix`, `array`, or `data.frame` corresponding to grid cell draws
  that will be aggregated to polygons.

- aggregation_table:

  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  Aggregation table linking pixels to polygon identifiers, created using
  [`build_aggregation_table()`](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md)

- aggregation_cols:

  (character vector, default 'polygon_id') Polygon identifiers to use
  for aggregation.

- method:

  (character, default 'mean') Aggregation method: one of 'mean', 'sum',
  'weighted.mean', or 'weighted.sum'.

- z_dimension:

  (vector, default NULL) If passing a `draws_matrix` with multiple sets
  of estimates, how should each layer be identified?

- z_dimension_name:

  (default 'z') The field name for the "z" dimension corresponding to
  each set of estimates contained in `draws_matrix`.

- weighting_raster:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  default NULL) The relative weighting of each whole pixel to the
  overall polygon value, for example, if calculating a
  population-weighted mean.

- na.rm:

  (bool, default TRUE) How to handle NA values in `draws_matrix` and
  `weighting_raster`.

## Value

Errors if checks fail; silently passes if checks pass

## See also

aggregate_draws_to_polygons
