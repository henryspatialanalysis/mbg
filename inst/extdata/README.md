# MBG data objects

The template files in this folder can be edited and used alongside the scripts in `inst/scripts`.

## Configuration file

The configuration file, `config.yaml`, includes settings related to the execution of the MBG scripts. The file is intended to be used with the `versioning` package, and more details about its structure are available in the [documentation](https://github.com/henryspatialanalysis/versioning) for that package.

After loading a config file, settings can be retrieved using the `$get()` method, directories can be retrieved using the `$get_dir_path()` method, and files can be read or written using the `$read()` and `$write()` methods. For example:

```
# Load config from file
path_to_config_file <- '/path/to/config.yaml'
config <- versioning::Config$new(path_to_config_file)

# Retrieve a top-level setting
config$get('indicator')

# Retrieve a nested setting
config$get('inla_settings', 'priors', 'range')

# Read the "adm2" file from the "shps" directory
adm_boundaries <- config$read("shps", "adm2")
# Save an object to the "adm_boundaries" file within the versioned "results" directory
config$write(adm_boundaries, "results", "adm_boundaries")
```


### Configuration file: settings

This section highlights configuration items that change model settings.

- `indicator` (string): Model indicator to run. This affects where the model execution script searches for the input data file and the title of certain diagnostic plots.
- `country` (string): Full name of the country where the model will be run. At the moment, the `mbg` package only supports single-country models. Should match the `ADM0_NAME` field in the modeling shapefile
- `iso3` (string): Capitalized ISO-3 code corresponding for `country`
- `year` (integer): Year when the model will be run. This affects how the input data is subset and which time-varying covariate files are used. At the moment, the `mbg` package only supports single-year models.
- `modeling_crs` (string): The coordinate reference system used for the input data, shapefile, and covariates. Typically `"EPSG:4326"`, the code for WGS84 unprojected lat-long.
- `draws_link_function`: Link function to connect model predictions (often estimated in a transformed space) with measurement space. Should be the name of an R function that matches the link function in `inla_settings/link`. For binomial logit-linked models, this should be set to `plogis` (the inverse logit function); for Poisson log-linked models, this should be set to `exp` (exponentiate, the inverse log function).
- `n_samples` (integer): Number of predictive draws to take at each pixel. Common settings are `100` or `250`
- `ui_width` (float < 1): Width of the uncertainty interval. Most commonly `0.95` for the 95% uncertainty interval, meaning that the "lower" estimates will be taken from the 2.5th percentile of draws, and "upper" estimates will be taken from the 97.5th percentile.
- `save_full_model` (logical): Should all model results be saved? Typically `TRUE` unless this is a test model or one of many out-of-sample runs.
- `shapefile_settings`:
  - `modeling_level` (int > 0): Most granular administrative level that the model will be aggregated to. Should be at least `1` (state/province level), and typically set to `2` (district/county level). The `ids` sub-list should contain administrative identifiers for all admin levels up to the modeling level---for example, if set to `2`, then `ids` should contain `adm0`, `adm1`, and `adm2` identifiers.
  - `ids`:
    - `adm0` (array of strings): Unique admin0 (country) identifiers. For example: `['ADM0_NAME', 'ADM0_ID']`
    - `adm1` (array of strings): Unique admin1 (state/province) identifiers. For example: `['ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE']`
    - `adm2` (array of strings): Unique admin2 (district/county) identifiers. For example: `['ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE', 'ADM2_NAME', 'ADM2_CODE']`
- `covariate_settings`:
  - `file_format` (string): File extension for covariate raster files. Often `'tif'` 
  - `add_intercept` (logical): Should an intercept covaraite containing all 1s be added to the list of loaded raster covariates?
- `pop_covariate_settings`:
  - `covariate` (string): Name of the population covariate to be loaded
  - `annual` (logical): Is the population covariate time-varying? Typically `TRUE`
  - `transform` (string): Name of the transformation to apply to the population covariate before returning. Typically `"identity"` (no transformation applied).
  - `normalize` (logical): Should the population surface be rescaled to N(0, 1) before use? Typically `FALSE` 
- `inla_settings`:
  - `family` (string): Model family to be used---for example, `"binomial"` or `"link"`. For more information, see the `INLA::INLA()` function [documentation](https://rdrr.io/github/andrewzm/INLA/man/inla.html).
  - `link` (string): Link function corresponding to `family`. Should also be reflected by `draws_link_function`.
  - `sample_size_field`: Field listing the total sample size for each data observation.
  - `priors`: Penalized complexity priors on the covariate fixed effects, the range of the Gaussian process, and the variance (standard deviation) of the Gaussian process.
    - `fixed_effects`: Prior on each fixed effect coefficient, excluding the intercept. Interpreted as the user-defined prior probability that each fixed effect coefficient will EXCEED a user-defined threshold.
      - `threshold` (numeric): UPPER threshold for the fixed effect coefficient to (possibly) exceed, for example, `3.0`
      - `prob_above` (numeric < 1): (Relatively small) probability that each fixed effect coefficient will exceed the threshold. Often set to `0.05` or `0.01`
    - `range`: Prior on the spatial range of the Gaussian process. Interpreted as the user-defined probability that the range will be BELOW a user-defined threshold.
      - `threshold` (numeric): LOWER threshold for the range to (possibly) go below. Set as a proportion of the largest spatial dimension across the modeling region: for example, if the study area is 500km across at its widest, a threshold of `0.2` is equivalent to a threshold of 100km
      - `prob_below` (numeric < 1): (Relatively small) probability that the range will go below the threshold. Often set to `0.05` or `0.01`
    - `sigma` (numeric): Prior on the standard deviation of the Gaussian process. Interpreted as the user-defined prior probability that sqrt(variance) of the Gaussian process will EXCEED a user-defined threshold. 
      - `threshold` (numeric): UPPER threshold for the GP standard deviation to (possibly) exceed.
      - `prob_above` (numeric < 1): (Relatively small) probability that the standard deviation will be above the threshold. Often set to `0.05` or `0.01`


### Configuration file: directory and file paths

The `directories` and `versions` sub-lists are special items that facilitate fast and versioned file reading/writing. Each item in the `directories` sub-list should be set up as follows:

```
directories:
    (directory name):
        path: "/path/to/directory"
        versioned: (TRUE or FALSE)
        files:
            (file name 1): "file_path_within_directory"
            (file name 2): "file_path_within_directory"
            ...
```

The `versions` sub-list should have entries corresponding to each directory where `versioned = TRUE`:

```
versions:
    (versioned directory 1): (version, e.g. "v1")
    (versioned directory 2): (version, e.g. "2023_10_26")
    ...
```

These sub-lists can be accessed using `config$get_dir_path()`, `config$read()`, and `config$write()`. For example, given the config:
```
directories:
    inputs:
        path: "/path/to/inputs"
        versioned: FALSE
        files:
            in_data: "a.csv"
    outputs:
        path: "/path/to/outputs"
        versioned: TRUE
        files:
            out_data: "b.csv"
versions:
    outputs: "v5"
```

Config methods can be used to easily load and save the files:

```
config$get_dir_path("inputs") # Not versioned - returns "/path/to/inputs"
config$get_dir_path("outputs) # Versioned - returns "/path/to/outputs/v5"

# Load "/path/to/inputs/a.csv" as a data.table
my_data <- config$read("inputs", "in_data") 

# Save an object to "/path/to/outputs/v5/b.csv"
config$write(my_data, "outputs", 'out_data")
```


## Example covariates

The covariates table includes information about where and how to load raster covariates. The fields include:

- `covariate` (string): Name of the raster covariate
- `include` (logical): Should the covariate be loaded in this model run? Allows for the user to easily include/exclude candidate covariates by toggling this field
-  `annual` (logical): Does the covariate take a different value in each year? Changes how the covariate loading function looks for the correct input file.
- `transform` (string): A transformation that should be applied to the covariate values before using them in regression. Can take the name of any common R function, and common candidates are `"identity"` (no transformation), `"log10"`, `"log1p"`, and `"sqrt"`.
- `normalize` (logical): Should the covariate be rescaled to approximately N(0, 1) before use? If `TRUE`, this rescaling happens after any other transformation.
