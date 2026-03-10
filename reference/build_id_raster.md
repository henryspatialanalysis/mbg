# Build ID raster

Build an ID raster matching the extent of a vector dataset

## Usage

``` r
build_id_raster(polygons, template_raster = NULL)
```

## Arguments

- polygons:

  [terra::SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  object. The polygons to be aggregated to

- template_raster:

  (optional)
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object. The template raster should contain and have the same CRS as
  the polygons. If template raster is `NULL`, the default, uses the
  default world template raster from
  [`make_world_template_raster()`](https://henryspatialanalysis.github.io/mbg/reference/make_world_template_raster.md).

## Value

ID raster. A
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object that minimally encloses the polygons

## Details

The ID raster will be used to build the aggregation table. Each pixel
has a unique integer value from 1 to the number of pixels in the ID
raster.

## Examples

``` r
if (FALSE) { # \dontrun{
  polygons <- sf::st_read(system.file('extdata/Benin_communes.gpkg', package = 'mbg'))
  build_id_raster(polygons)
} # }
```
