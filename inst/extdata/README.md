# MBG data objects

The template files in this folder can be edited and used alongside the scripts in `inst/scripts`.

## Configuration file

The configuration file, `config.yaml`, includes settings related to the execution of the MBG scripts. The file is intended to be used with the `versioning` package, and more details about its structure are available in the [documentation](https://github.com/henryspatialanalysis/versioning) for that package.

## Example covariates

The covariates table includes information about where and how to load raster covariates. The fields include:

- `covariate` (string): Name of the raster covariate
- `include` (logical): Should the covariate be loaded in this model run? Allows for the user to easily include/exclude candidate covariates by toggling this field
-  `annual` (logical): Does the covariate take a different value in each year? Changes how the covariate loading function looks for the correct input file.
- `transform` (string): A transformation that should be applied to the covariate values before using them in regression. Can take the name of any common R function, and common candidates are `"identity"` (no transformation), `"log10"`, `"log1p"`, and `"sqrt"`.
- `normalize` (logical): Should the covariate be rescaled to approximately N(0, 1) before use? If `TRUE`, this rescaling happens after any other transformation.
