# Aggregate a raster to polygons: validation

Data validation for aggregate_raster_to_polygons

## Usage

``` r
aggregate_raster_to_polygons_validation(
  data_raster,
  aggregation_table,
  aggregation_cols,
  method,
  aggregated_field,
  z_dimension,
  z_dimension_name,
  weighting_raster,
  na.rm
)
```

## Arguments

- data_raster:

  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  containing data to be aggregated to polygons.

- aggregation_table:

  Aggregation table linking pixels to polygon identifiers, created using
  [`build_aggregation_table()`](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md)

- aggregation_cols:

  (character vector, default 'polygon_id') Polygon identifiers to use
  for aggregation.

- method:

  (character, default 'mean') Aggregation method: one of 'mean', 'sum',
  'weighted.mean', or 'weighted.sum'.

- aggregated_field:

  (character, default 'data') Name of the aggregated raster values in
  the output table.

- z_dimension:

  (vector, default NULL) If passing a `data_raster` with multiple
  layers, how should each layer be identified?

- z_dimension_name:

  (default 'z') The field name for the "z" dimension corresponding to
  each layer of the `data_raster`.

- weighting_raster:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  default NULL) The relative weighting of each whole pixel to the
  overall polygon value, for example, if calculating a
  population-weighted mean.

- na.rm:

  (bool, default TRUE) How to handle NA pixels in `data_raster` and
  `weighting_raster`.

## Value

Errors if checks fail; silently passes if checks pass

## See also

aggregate_raster_to_polygons
