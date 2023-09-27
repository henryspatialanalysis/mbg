## #######################################################################################
##
## MBG package workflow
##
## AUTHOR: Nathaniel Henry, nat@henryspatialanalysis.com
## CREATED: 4 September 2023
## PURPOSE: MBG workflow using packages `mbg`, `versioning`, and `pixel2poly`
##
## #######################################################################################

## 00) SETTINGS ------------------------------------------------------------------------->

# Default locations of repositories and configuration file, if not specified below

DEFAULT_REPOS_PATH <- '~/repos'
DEFAULT_CONFIG_PATH <- file.path(DEFAULT_REPOS_PATH, 'mbg/inst/extdata/example_config.yaml')


# You can optionally update certain important settings by run without editing the config.
#  All of these have default values that will be used if they are not set. These settings
#  include:
#   - `repos_path`: Path to the folder that contains all repositores. Defaults to "~/repos"
#   - `config_path`: Path to the config.yaml file to use for this model run. Defaults to
#       "<repos_path>/mbg/inst/extdata/example_config.yaml"
#   - `indicator`: Model indicator to run. Defaults to config value.
#   - `country`: ISO3 code to run. Defaults to config value.
#   - `year`: Model year. Defaults to config value.
#   - `results_version`: Determines the time-stamped folder where results will be saved.
#       Defaults to config value.
#
# If you are running interactively, you can optionally update these settings on lines 
#  45-51. Except for `repos_path` and `config_path`, which always need to be set, leaving
#  any of these arguments blank will default to the settings already listed in the config.
#
# If you are running on the command line, you can pass these settings using named command
#  line arguments. If an argument is not passed, the default settings will be used 
#  instead. For example, this call to the script will update the "indicator" and "year",
#  but use default settings otherwise:
#
# Rscript --vanilla mbg_package_workflow.R --indicator modern_cpr --year 2021

if(interactive()){
  # If running the script interactively, you can optionally update settings
  run_specific_settings <- list(
    repos_path = DEFAULT_REPOS_PATH,
    config_path = DEFAULT_CONFIG_PATH,
    indicator = 'stunted',
    iso3 = 'TZA',
    country = 'Tanzania',
    year = 2022,
    results_version = '20230927'
  )
} else {
  library(argparse)
  parser <- argparse::ArgumentParser()
  parser$add_argument("--repos_path", type = 'character', default = DEFAULT_REPOS_PATH)
  parser$add_argument("--config_path", type = 'character', default = DEFAULT_CONFIG_PATH)
  parser$add_argument("--indicator", type = 'character', default = NULL)
  parser$add_argument("--iso3", type = 'character', default = NULL)
  parser$add_argument("--country", type = 'character', default = NULL)
  parser$add_argument("--year", type = 'integer', default = NULL)
  parser$add_argument("--results_version", type = "character", default = NULL)
  run_specific_settings <- parser$parse_args(
    commandArgs(trailingOnly = TRUE)
  )
}


## 01) SETUP ---------------------------------------------------------------------------->

# Load standard packages
load_packages <- c(
  'assertthat', 'data.table', 'glue', 'INLA', 'rgeoboundaries', 'sf', 'terra', 'tictoc'
)
invisible(lapply(load_packages, library, character.only = TRUE))
tictoc::tic("Full script execution")

# Load custom packages from subfolders within the `repos_path` directory
assertthat::assert_that(dir.exists(run_specific_settings$repos_path))
for(custom_package in c('versioning', 'pixel2poly', 'mbg')){
  custom_package_dir <- file.path(run_specific_settings$repos_path, custom_package)
  assertthat::assert_that(dir.exists(custom_package_dir))
  devtools::load_all(custom_package_dir)
}

# Set the custom results version, if one was passed 
if(is.null(run_specific_settings$results_version)){
  custom_versions <- NULL
} else {
  custom_versions <- list(results = run_specific_settings$results_version)
}

# Load the configuration object
config <- versioning::Config$new(
  config_list = run_specific_settings$config_path,
  versions = custom_versions
)

# Update the config with run-specific settings, if they were passed
for(setting_name in c('indicator', 'iso3', 'country', 'year')){
  custom_value <- run_specific_settings[[setting_name]]
  if(!is.null(custom_value)){
    config$config_list[[setting_name]] <- custom_value
  }
}
# The name of the input file is "<raw_data directory>/<indicator>.csv"
config$config_list$directories$raw_data$files <- list(
  input_data = paste0(config$get('indicator'), '.csv')
)

# Create the results directory
results_dir <- config$get_dir_path('results')
dir.create(results_dir, recursive = T, showWarnings = F)

# Save a copy of the config in the results directory
config$write_self('results')


## 02) PREPARE INPUT DATA --------------------------------------------------------------->

# Load country boundaries
country <- config$get('country')
adm_level <- paste0('adm', config$get('shapefile_settings', 'modeling_level'))
shp_id_cols <- config$get('shapefile_settings', 'ids', adm_level)
shp_query <- glue::glue(
  'SELECT {paste(shp_id_cols, collapse = ",")} ',
  'FROM {tools::file_path_sans_ext(config$get("directories", "shps", "files", adm_level))} ',
  'WHERE ADM0_NAME = \'{config$get("country")}\''
)
adm_boundaries <- config$read(dir_name = "shps", file_name = adm_level, query = shp_query)
adm_boundaries$polygon_id <- seq_len(nrow(adm_boundaries))

# Create the ID raster from the admin2 spatial object
id_raster <- pixel2poly::build_id_raster(polygons = terra::vect(adm_boundaries))

# Load a list of covariates
covariates_list <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('covariates'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariate_settings', 'file_format'),
  add_intercept = config$get('covariate_settings', 'add_intercept')
)
# Load the population raster
# Used for population-weighted aggregation from grid cell to admin
population_raster <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('pop_covariate'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariate_settings', 'file_format'),
  add_intercept = FALSE
)[[1]]

# Load input data from the `raw_data` directory
# This file path was set based on the indicator on line 105
input_data <- (
  config$read("raw_data", "input_data")
  [(year == config$get("year")) & (country == config$get("iso3")), ] |>
  setnames(
    old = c('longitude', 'latitude', 'N', config$get('indicator')),
    new = c('x', 'y', 'samplesize', 'indicator'),
    skip_absent = TRUE
  )
)
if(nrow(input_data) == 0) stop(
  "After subsetting to data from ", config$get("country"), " in ", config$get("year"),
  ", no rows of data remain."
)


## 03) PREPARE AND RUN INLA MODEL ------------------------------------------------------->

inla_inputs_list <- mbg::prepare_inla_data_stack(
  input_data = input_data,
  id_raster = id_raster,
  covariates = covariates_list
)

inla_fitted_model <- mbg::fit_inla_model(
  formula = inla_inputs_list$formula_string,
  data_stack = inla_inputs_list$inla_data_stack,
  spde = inla_inputs_list$spde,
  family = config$get('inla_settings', 'family'),
  link = config$get('inla_settings', 'link')
)


## 05) GENERATE PREDICTIONS ------------------------------------------------------------->

# Predictions by grid cell
grid_cell_predictions <- mbg::generate_cell_draws_and_summarize(
  inla_model = inla_fitted_model,
  inla_mesh = inla_inputs_list$mesh,
  n_samples = config$get("n_samples"),
  id_raster = id_raster,
  covariates = covariates_list,
  inverse_link_function = config$get('draws_link_function'),
  ui_width = config$get('ui_width')
)

# Prepare objects to store admin draws and summaries
max_adm_level <- config$get('shapefile_settings', 'modeling_level')
max_adm_level_label <- paste0('adm', max_adm_level)
all_adm_levels <- seq(0, max_adm_level)
adm_draws_list <- adm_summaries_list <- vector('list', length = length(all_adm_levels))
names(adm_draws_list) <- names(adm_summaries_list) <- paste0('adm', all_adm_levels)
draw_fields <- paste0('draw_', seq_len(config$get('n_samples')))

# Aggregate to the most detailed admin units
message("Aggregating draws at the ", max_adm_level_label, " level.")
aggregation_table <- pixel2poly::build_aggregation_table(
  polygons = terra::vect(adm_boundaries),
  id_raster = id_raster,
  polygon_id_field = 'polygon_id'
)
agg_cols <- config$get("shapefile_settings", "ids", max_adm_level_label)
detailed_adm_draws <- pixel2poly::aggregate_draws_to_polygons(
  draws_matrix = grid_cell_predictions$cell_draws,
  aggregation_table = aggregation_table,
  aggregation_cols = agg_cols,
  method = 'weighted.mean',
  weighting_raster = population_raster
)
admin_pop <- pixel2poly::aggregate_raster_to_polygons(
  data_raster = population_raster,
  aggregation_table = aggregation_table,
  aggregation_cols = agg_cols,
  method = 'sum',
  aggregated_field = 'population'
)
detailed_adm_draws[admin_pop, population := i.population, on = agg_cols]
adm_draws_list[[max_adm_level_label]] <- data.table::copy(detailed_adm_draws)

# Aggregate to less-detailed admin units
for(higher_level in setdiff(names(adm_draws_list), max_adm_level_label)){
  message("Aggregating draws at the ", higher_level, " level.")
  agg_cols <- config$get("shapefile_settings", "ids", higher_level)
  adm_draws_list[[higher_level]] <- detailed_adm_draws[
    , c(
      list(population = sum(population)),
      lapply(.SD, weighted.mean, w = population)
    ),
    .SDcols = draw_fields,
    by = agg_cols
  ]
}

# Summarize admin draws at all levels
for(adm_level in names(adm_summaries_list)){
  id_fields <- config$get("shapefile_settings", "ids", adm_level)
  adm_summary <- mbg::summarize_draws(
    draws = adm_draws_list[[adm_level]],
    id_fields = id_fields,
    draw_fields = draw_fields,
    ui_width = config$get("ui_width")
  )
  adm_summary[adm_draws_list[[adm_level]], population := i.population, on = id_fields]
  adm_summaries_list[[adm_level]] <- adm_summary
}


## 06) SAVE RESULTS --------------------------------------------------------------------->

config$write(adm_boundaries, 'results', 'adm_boundaries')
config$write(id_raster, 'results', 'id_raster')
config$write(terra::rast(covariates_list), 'results', 'covariate_rasters')
config$write(input_data, 'results', 'formatted_input_data')
config$write(inla_inputs_list, 'results', 'inla_data_stack')
if(config$get('save_full_model')) config$write(inla_fitted_model, 'results', 'inla_model')
config$write(grid_cell_predictions$parameter_draws, 'results', 'parameter_draws')
config$write(grid_cell_predictions$cell_draws, 'results', 'cell_draws')
config$write(grid_cell_predictions$cell_pred_mean, 'results', 'cell_pred_mean')
config$write(grid_cell_predictions$cell_pred_lower, 'results', 'cell_pred_lower')
config$write(grid_cell_predictions$cell_pred_upper, 'results', 'cell_pred_upper')
for(adm_level in names(adm_draws_list)){
  config$write(adm_draws_list[[adm_level]], 'results', paste0(adm_level, '_draws'))
  config$write(adm_summaries_list[[adm_level]], 'results', paste0(adm_level, '_summary_table'))
}

# End script timer
tictoc::toc()
