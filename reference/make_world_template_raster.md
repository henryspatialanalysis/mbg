# Make world template raster

Create a template raster for the world with approximately 5x5km
resolution at the equator, matching many common raster covariates for
health.

## Usage

``` r
make_world_template_raster()
```

## Value

[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object matching the specifications above

## Details

The raster has the following specifications:

- 4320 rows, 8640 columns

- Resolution: 0.04166667 decimal degrees, approx. 5km at the equator

- CRS: WGS 1984 unprojected latitude/longitude,
  `terra::crs('EPSG:4326')`

- Values: All NA. Used exclusively for creating a shapefile-specific ID
  raster
