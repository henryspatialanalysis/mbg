# Insert values into a raster

Insert a vector or matrix of values into an ID spatRaster

## Usage

``` r
values_to_raster(x, id_raster)
```

## Arguments

- x:

  Vector, matrix, data.frame, or data.table of values that will be
  inserted into the ID raster. The length of x must be exactly divisible
  by `sum(!is.na(terra::values(id_raster)))`. Data.frames are converted
  to matrices, and then matrices are converted to vectors using
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) and
  [`as.vector()`](https://rdrr.io/r/base/vector.html) respectively
  before processing. For that reason, data.frames should only contain
  fields with values to be inserted (such as a data.frame of draws).

- id_raster:

  ID raster showing the outline of the study area, created using
  [`build_id_raster()`](https://henryspatialanalysis.github.io/mbg/reference/build_id_raster.md).
  Should have 1 layer.

## Value

SpatRaster with the same outline as the ID raster and (# values / \#
non-NA pixels in the ID raster) layers.

## Details

The length of the vector or matrix must be a multiple of the number of
non-NA pixels in the ID raster. Values from the vector/matrix are then
inserted into the non-NA pixels of the spatRaster.

## See also

[`build_id_raster()`](https://henryspatialanalysis.github.io/mbg/reference/build_id_raster.md)

## Examples

``` r
# Example ID raster with 10 rows and 10 columns, and 99 valid pixels
example_id_raster <- terra::rast(matrix(c(seq_len(99), NA), nrow = 10))
# Inserting 99 values yields a spatRaster with 1 layer
mbg::values_to_raster(stats::rnorm(99), example_id_raster)
#> class       : SpatRaster 
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> name        :     lyr.1 
#> min value   : -2.938978 
#> max value   :  2.101253 
# Inserting 99 * 3 values yields a spatRaster with 3 layers
mbg::values_to_raster(seq_len(99 * 3), example_id_raster)
#> class       : SpatRaster 
#> size        : 10, 10, 3  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
#> source(s)   : memory
#> names       : lyr.1, lyr.1, lyr.1 
#> min values  :     1,   100,   199 
#> max values  :    99,   198,   297 
# Trying to insert values with length not divisible by 99 yields an error
try(mbg::values_to_raster(seq_len(100), example_id_raster))
#> Error : The length of x must be exactly divisible by the number of non-NA pixels in the ID raster.
```
