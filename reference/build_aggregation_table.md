# Build aggregation table

Build a table to quickly aggregate from pixels to polygons

## Usage

``` r
build_aggregation_table(polygons, id_raster, polygon_id_field, verbose = FALSE)
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

  (`character(1)`) Unique identifier field in `polygons`.

- verbose:

  (`logical(1)`, default FALSE) Show progress for building aggregation
  rows for each polygon?

## Value

data.table with fields:

- polygon_id: Unique polygon identifier

- pixel_id: unique pixel ID from the ID raster

- masked_pixel_id: Index counting only non-NA pixels from the ID raster

- area_fraction: fraction of the pixel area falling within this polygon

- Merged fields from the table of polygons

## See also

calculate_pixel_fractions_single_polygon()

## Examples

``` r
if (FALSE) { # \dontrun{
  polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
  id_raster <- build_id_raster(polygons)
  aggregation_table <- build_aggregation_table(
    polygons, id_raster, polygon_id_field = 'commune_code'
  )
} # }
```
