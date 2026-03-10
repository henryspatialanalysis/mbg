# Load covariates

Load covariates

## Usage

``` r
load_covariates(
  directory,
  covariates_table,
  id_raster,
  year = NULL,
  file_format = "tif",
  add_intercept = FALSE,
  check_previous_years = 10
)
```

## Arguments

- directory:

  Directory containing all covariate sub-directories

- covariates_table:

  `data.frame` containing at least the following fields:

  - 'covariate': (character): Name of the covariate

  - 'annual': (logical) Does the covariate vary by year? If so, look for
    the `year` in the name of the file.

  - 'transform': (character) Name of a function to use to transform the
    covariate. Common options include 'identity' (no transformation),
    'sqrt', 'abs', and 'log1p'

  - 'normalize': (logical) Should the covariate be rescaled to have mean
    0 and standard deviation 1 across all pixels in the study area?
    Generally should be set to TRUE for predictive covariates.

- id_raster:

  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with non-NA pixels delineating the extent of the study area

- year:

  (`numeric`, default NULL) Year of data to for time-varying covariates.
  If NULL, the default, uses the current year.

- file_format:

  (`character`, default 'tif') File format for the raster covariate
  data. Used to search for the input file within the proper containing
  folder.

- add_intercept:

  (`logical`, default FALSE) Should a covariate called "intercept", a
  raster object with 1s in all required cells, be placed at the start of
  the returned covariates list?

- check_previous_years:

  (`integer` \> 0, default 10) If annual data is not found in this year,
  how many previous years should be checked? If 0, will not check any
  previous years.

## Value

A named list of formatted covariates. Each list item is a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with one layer and the same dimensions as the `id_raster`

## Details

Load a list of covariates from a specified directory structure
