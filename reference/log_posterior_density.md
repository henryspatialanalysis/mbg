# Generate log posterior predictive density from a geostatistical surface onto point data

Generate log posterior predictive density from a geostatistical surface
onto point data

## Usage

``` r
log_posterior_density(draws, validation_data, id_raster, na.rm = FALSE)
```

## Arguments

- draws:

  (`matrix`) A predictive draw matrix, where each row corresponds to a
  pixel in the `id_raster` and each column corresponds to one sampled
  estimate of the outcome.

- validation_data:

  (`data.frame`) Table containing at least the following fields:  

  - x (`numeric`) location x position, in the same projection as
    `id_raster`  

  - y (`numeric`) location y position, in the same projection as
    `id_raster`  

  - indicator (`integer`) The number of events in the population  

  - samplesize (`integer`) The total population, denominator for
    `indicator`

- id_raster:

  ([terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html))
  Raster showing the sample study area, created using
  [build_id_raster](https://henryspatialanalysis.github.io/mbg/reference/build_id_raster.md).

- na.rm:

  (`logical(1)`, default FALSE) Should NA values be omitted from the LPD
  calculation?

## Value

(`numeric(1)`) Log predictive density of the validation data given the
draw estimates.

## Details

Calculated across draws. Requires an ID raster to match each point
observation to a set of draws. Assumes binomial data.

For examples, see
[`vignette('model-comparison', package = 'mbg')`](https://henryspatialanalysis.github.io/mbg/articles/model-comparison.md)
