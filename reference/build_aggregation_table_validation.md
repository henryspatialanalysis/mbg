# Validation: Build aggregation table

Input data validation for
[`build_aggregation_table()`](https://henryspatialanalysis.github.io/mbg/reference/build_aggregation_table.md)

## Usage

``` r
build_aggregation_table_validation(
  polygons,
  id_raster,
  polygon_id_field,
  polygon_ids
)
```

## Arguments

- polygons:

  [terra::SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  object. Should contain a unique ID field.

- id_raster:

  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object. ID raster created by
  [`build_id_raster()`](https://henryspatialanalysis.github.io/mbg/reference/build_id_raster.md)
  for the polygons object. Should have the same CRS as `polygons` and
  completely cover it.

- polygon_id_field:

  (`character`) Unique identifier field in `polygons`.

- polygon_ids:

  (vector, any type) Polygon identifiers from `polygon_id_field`.

## Value

Errors if checks fail; silently passes if checks pass

## See also

build_aggregation_table
