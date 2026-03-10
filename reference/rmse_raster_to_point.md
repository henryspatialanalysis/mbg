# Generate RMSE from an estimated raster surface and point data

Generate RMSE from an estimated raster surface and point data

## Usage

``` r
rmse_raster_to_point(estimates, validation_data, outcome_field, na.rm = FALSE)
```

## Arguments

- estimates:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))
  Raster surface containing point estimates. This could also be the mean
  surface of a Bayesian geostatistical model

- validation_data:

  (`data.frame`)  
  Table containing at least the following fields:  

  - x (`numeric`) location x position, in the same projection as
    `estimates`  

  - y (`numeric`) location y position, in the same projection as
    `estimates`  

  - (Outcome field) See below

- outcome_field:

  (`character(1)`) Column in `validation_data` containing the values
  that should be compared against the `estimates` raster surface.

- na.rm:

  (`logical(1)`, default FALSE) Should NA values be dropped from the
  RMSE calculation?

## Value

A single number giving RMSE between the point data and estimates raster.

## Details

For examples, see
[`vignette('model-comparison', package = 'mbg')`](https://henryspatialanalysis.github.io/mbg/articles/model-comparison.md)
