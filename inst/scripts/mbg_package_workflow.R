## #######################################################################################
##
## MBG package workflow
##
## AUTHOR: Nathaniel Henry, nat@henryspatialanalysis.com
## CREATED: 4 September 2023
## PURPOSE: Example MBG workflow using packages `mbg`, `versioning`, and `pixel2poly`
##
## #######################################################################################

## SETTINGS

# Eventually, only these settings will be changed for each script iteration. All other
#  settings can be changed in the config file.
custom_versions <- list(results = '20230904')


## 01) SETUP ---------------------------------------------------------------------------->

# Load standard packages
load_packages <- c('data.table', 'INLA', 'rgeoboundaries', 'sf', 'terra', 'tictoc')
invisible(lapply(load_packages, library, character.only = TRUE))

# Load custom packages
custom_pkgs <- c('versioning', 'pixel2poly', 'mbg')
load_locally <- TRUE
if(load_locally){
  repos_path <- '~/repos'
  lapply(file.path(repos_path, custom_pkgs), devtools::load_all) |> invisible()
} else {
  github_acct <- 'henryspatialanalysis'
  for(custom_pkg in custom_pkgs){
    remotes::install_github(file.path(github_acct, custom_pkg))
    library(custom_pkg, character.only = T)
  }
}

# Make config and create results directory
config_file_path <- system.file('extdata/example_config.yaml', package = 'mbg')

config <- versioning::Config$new(
  config_list = config_file_path,
  versions = custom_versions
)
results_dir <- config$get_dir_path('results')
dir.create(results_dir, recursive = T, showWarnings = F)
# Save a copy of the config in the results directory
config$write_self('results')


## 02) PREPARE INPUT DATA --------------------------------------------------------------->

# Load country boundaries
country <- config$get('country')
boundaries_list <- list(
  adm0 = rgeoboundaries::geoboundaries(country, adm_lvl = 0),
  adm1 = rgeoboundaries::geoboundaries(country, adm_lvl = 1),
  adm2 = rgeoboundaries::geoboundaries(country, adm_lvl = 2)
)
# Add a `polygon_id` to all boundaries objects
for(ii in seq_along(boundaries_list)){
  boundaries_list[[ii]]$polygon_id <- seq_len(nrow(boundaries_list[[ii]]))
}

# Create the ID raster from the admin2 spatial object
id_raster <- pixel2poly::build_id_raster(
  polygons = terra::vect(boundaries_list$adm2)
)

# Load a list of covariates
covariates_list <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('covariate_settings'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariates_file_format'),
  add_intercept = config$get('covariates_add_intercept')
)
# Load the population raster
# Used for population-weighted aggregation from grid cell to admin
population_raster <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('pop_covariate_settings'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariates_file_format'),
  add_intercept = FALSE
)[[1]]

# Load input data
input_data <- (
  data.table::as.data.table(get(data(gambia, package = 'geoR')))
  [, .(indicator = sum(pos), samplesize = .N), by = .(x, y)]
)
# Convert to working CRS (typically unprojected lat-long)
latlong <- (
  sf::st_as_sf(input_data, coords = c('x', 'y'), crs = sf::st_crs('+proj=utm +zone=28')) |>
  sf::st_transform(crs = sf::st_crs(config$get('crs'))) |>
  sf::st_coordinates() |>
  as.data.table()
)
input_data$x <- latlong$X
input_data$y <- latlong$Y


## 03) PREPARE AND RUN INLA MODEL ------------------------------------------------------->

inla_inputs_list <- mbg::prepare_inla_data_stack(
  input_data = input_data,
  id_raster = id_raster,
  covariates = covariates_list
)

inla_fitted_model <- mbg::fit_inla_model(
  formula = as.formula(inla_inputs_list$formula_string),
  data_stack = inla_inputs_list$inla_data_stack,
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

# Aggregate to predictions by admin unit
admin_levels <- names(boundaries_list)
admin_predictions_list <- vector('list', length = length(admin_levels))
names(admin_predictions_list) <- admin_levels
for(admin_level in admin_levels){
  aggregation_table <- pixel2poly::build_aggregation_table(
    polygons = terra::vect(boundaries_list[[admin_level]]),
    id_raster = id_raster,
    polygon_id_field = 'polygon_id'
  )
  admin_draws <- pixel2poly::aggregate_draws_to_polygons(
    draws_matrix = grid_cell_predictions$cell_draws,
    aggregation_table = aggregation_table,
    aggregation_cols = config$get('admin_id_fields'),
    method = 'weighted.mean',
    weighting_raster = population_raster
  )
  admin_summaries <- mbg::summarize_draws(
    draws = admin_draws,
    id_fields = config$get('admin_id_fields'),
    ui_width = config$get('ui_width')
  )
  admin_predictions_list[[admin_level]] <- list(
    admin_draws = admin_draws,
    admin_summaries = admin_summaries
  )
}


## 06) SAVE RESULTS --------------------------------------------------------------------->

config$write(boundaries_list$adm0, 'results', 'adm0_shp')
config$write(boundaries_list$adm1, 'results', 'adm1_shp')
config$write(boundaries_list$adm2, 'results', 'adm2_shp')
config$write(id_raster, 'results', 'id_raster')
config$write(terra::rast(covariates_list), 'results', 'covariate_rasters')
config$write(inla_inputs_list, 'results', 'inla_data_stack')
config$write(inla_fitted_model, 'results', 'inla_model')
config$write(grid_cell_predictions$parameter_draws, 'results', 'parameter_draws')
config$write(grid_cell_predictions$cell_draws, 'results', 'cell_draws')
config$write(grid_cell_predictions$cell_pred_mean, 'results', 'cell_pred_mean')
config$write(grid_cell_predictions$cell_pred_lower, 'results', 'cell_pred_lower')
config$write(grid_cell_predictions$cell_pred_upper, 'results', 'cell_pred_upper')
config$write(admin_predictions_list$adm2$admin_draws, 'results', 'adm2_draws')
config$write(admin_predictions_list$adm2$admin_summaries, 'results', 'adm2_summary_table')
config$write(admin_predictions_list$adm1$admin_draws, 'results', 'adm1_draws')
config$write(admin_predictions_list$adm1$admin_summaries, 'results', 'adm1_summary_table')
config$write(admin_predictions_list$adm0$admin_draws, 'results', 'adm0_draws')
config$write(admin_predictions_list$adm0$admin_summaries, 'results', 'adm0_summary_table')
