## #######################################################################################
##
## Visualize MBG results
##
## AUTHOR: Nathaniel Henry, nat@henryspatialanalysis.com
## CREATED: 11 September 2023
## PURPOSE: Visualize outputs of mbg_package_workflow.R
##
## #######################################################################################

## 00) SETTINGS ------------------------------------------------------------------------->

REPOS_PATH_DEFAULT <- '~/repos'

# The script takes one required argument:
#  - `config_path`: the path to the `config.yaml` configuration file that was saved in
#    mbg_package_workflow.R
#
# The script also takes one optional argument, which has a default:
#  - `repos_path`: the path to the directory that contains the `versioning` custom package

if(interactive()){
  config_path <- '~/temp_data/geostats/mbg_results/20231019/config.yaml'
  repos_path <- REPOS_PATH_DEFAULT
} else {
  library(argparse)
  parser <- argparse::ArgumentParser()
  parser$add_argument("--config_path", type = 'character')
  parser$add_argument("--repos_path", type = 'character', default = REPOS_PATH_DEFAULT)
  # Parse command line arguments
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  config_path <- args$config_path
  repos_path <- args$repos_path
}


## 01) SETUP ---------------------------------------------------------------------------->

# Load standard packages
load_packages <- c(
  'data.table', 'ggplot2', 'glue', 'patchwork', 'RColorBrewer', 'sf', 'viridisLite', 'INLA'
)
invisible(lapply(load_packages, library, character.only = TRUE))

# Load the versioning package based on the `repos_path`
versioning_dir <- file.path(repos_path, 'versioning')
assertthat::assert_that(dir.exists(versioning_dir))
devtools::load_all(versioning_dir)

# Load the configuration object
config <- versioning::Config$new(config_list = config_path)

# Create the visualization directory as a sub-folder of the results directory
results_dir <- config$get_dir_path('results')
viz_dir <- file.path(results_dir, 'viz')
dir.create(viz_dir)

# Load data for visualization
# See the config.yaml file for the list of all results files that can be loaded
modeling_level <- config$get('shapefile_settings', 'modeling_level')
adm_boundaries <- config$read('results', 'adm_boundaries')
summary_table <- config$read('results', glue::glue('adm{modeling_level}_summary_table'))
mean_raster <- config$read('results', 'cell_pred_mean')
lower_raster <- config$read('results', 'cell_pred_lower')
upper_raster <- config$read('results', 'cell_pred_upper')
pop_raster <- config$read('results', 'pop_raster')
input_data <- config$read('results', 'formatted_input_data')
mesh <- config$read('results', 'inla_data_stack')$mesh

# Set some global variables that will be used in labels several times
indicator <- config$get('indicator')
country <- config$get('country')
year <- config$get('year')

# Set the color scheme for mean/lower/upper values
low_is_better <- TRUE
color_scheme <- RColorBrewer::brewer.pal(n = 9, name = 'Spectral')
if(low_is_better) color_scheme <- rev(color_scheme)

# Shortcut for element_blank()
eb <- ggplot2::element_blank()
map_theme <- ggplot2::theme_minimal() + 
  ggplot2::theme(axis.ticks = eb, axis.text = eb, panel.grid = eb)

## 02) Plot the mesh -------------------------------------------------------------------->

pdf(file.path(viz_dir, 'inla_mesh.pdf'), height = 10, width = 10)
plot(mesh)
dev.off()


## 03) Plot the summary rasters --------------------------------------------------------->

raster_to_table <- function(rr, type){
  r_table <- as.data.frame(rr, xy = T)
  names(r_table) <- c('x','y','value')
  r_table$type = type
  return(r_table)
}
raster_table_full <- rbind(
  raster_to_table(mean_raster, type = '1. Mean'),
  raster_to_table(lower_raster, type = '2. Lower'),
  raster_to_table(upper_raster, type = '3. Upper')
)
value_limits <- quantile(na.omit(raster_table_full$value), probs = c(0.05, 0.95))

raster_fig <- ggplot() + 
  facet_wrap('type', nrow = 1) +
  geom_raster(data = raster_table_full, aes(x = x, y = y, fill = value)) + 
  geom_sf(data = adm_boundaries, fill = NA, linewidth = 0.05, color = '#444444') + 
  scale_fill_gradientn(
    colors = color_scheme, limits = value_limits, oob = scales::squish,
    labels = scales::percent
  ) +
  labs(
    title = glue::glue("Summary rasters: {indicator} in {country}, {year}"),
    x = '', y = '', fill = indicator
  ) +
  map_theme
pdf(file.path(viz_dir, 'summary_rasters.pdf'), height = 6, width = 10)
print(raster_fig)
dev.off()


## 03) Plot the summary administrative data --------------------------------------------->

adm_ids <- config$get('shapefile_settings', 'ids', paste0('adm', modeling_level))
adm_data_for_plotting <- merge(
  x = adm_boundaries,
  y = melt(summary_table, id.vars = adm_ids),
  by = adm_ids
) |> merge(
  y = data.frame(
    variable = c('mean','lower','upper'),
    type = c('1. Mean', '2. Lower', '3. Upper')
  ),
  by = 'variable'
)

adm_fig <- ggplot() + 
  facet_wrap('type', nrow = 1) +
  geom_sf(
    data = adm_data_for_plotting,
    aes(fill = value), linewidth = 0.25, color = '#444444'
  ) + 
  scale_fill_gradientn(colors = color_scheme, labels = scales::percent) +
  labs(
    title = glue::glue("Summary estimates: {indicator} in {country}, {year}"),
    x = '', y = '', fill = indicator
  ) +
  map_theme
pdf(file.path(viz_dir, 'summary_admin_estimates.pdf'), height = 6, width = 10)
print(adm_fig)
dev.off()


## 04) Plot the input point data -------------------------------------------------------->

input_data[, rate := indicator / samplesize ]
point_color_limits <- quantile(na.omit(input_data$rate), probs = c(0.05, 0.95))

data_point_fig <- ggplot() + 
  geom_sf(data = adm_boundaries, fill = NA, linewidth = 0.25, color = "#444444") + 
  geom_point(
    data = input_data,
    aes(x = x, y = y, fill = rate, size = samplesize),
    shape = 21, alpha = .8
  ) + 
  scale_fill_gradientn(
    colors = color_scheme, labels = scales::percent, limits = point_color_limits,
    oob = scales::squish
  ) + 
  labs(
    title = glue::glue("Raw data: {indicator} in {country}, {year}"),
    x = '', y = '', fill = stringr::str_to_title(indicator), size = 'Number\nsampled'
  ) +
  map_theme

pdf(file.path(viz_dir, 'input_data_map.pdf'), height = 8, width = 8)
plot(data_point_fig)
dev.off()
png(file.path(viz_dir, 'input_data_map.png'), height = 1600, width = 1600, res = 200)
plot(data_point_fig)
dev.off()

## 04) Plot mean vs. uncertainty bivariate maps ----------------------------------------->

bv_color_df <- data.table::CJ(uncertainty_group = 0:2, value_group = 0:2)
bv_color_df[, idx := as.character(.I) ]
bv_color_vec <- c(
  "#e8e8e8", "#e4acac", "#c85a5a", "#b0d5df", "#ad9ea5", "#985356",
  "#64acbe", "#627f8c", "#574249"
)
names(bv_color_vec) <- bv_color_df$idx
color_scheme_fig <- ggplot() +
  geom_raster(data = bv_color_df, aes(x = value_group, y = uncertainty_group, fill = idx)) + 
  labs(title = '', x = 'Mean estimates', y = 'Uncertainty') +
  scale_fill_manual(values = bv_color_vec) +
  guides(fill = 'none') +
  theme_minimal() + 
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
    axis.title = element_text(size = 12, face = 'italic'),
    plot.background = element_rect(fill='transparent', color=NA)
  )

bv_keep_fields <- c(adm_ids, 'idx')
bv_data_dt <- (copy(summary_table)
  [, ui := upper - lower ]
  [, value_group := 2 ]
  [ mean < quantile(mean, 2/3, na.rm= T), value_group := 1]
  [ mean < quantile(mean, 1/3, na.rm= T), value_group := 0 ]
  [, uncertainty_group := 2 ]
  [ ui < quantile(ui, 2/3, na.rm= T), uncertainty_group := 1 ]
  [ ui < quantile(ui, 1/3, na.rm= T), uncertainty_group := 0 ]
  [bv_color_df, idx := i.idx, on = c('value_group', 'uncertainty_group')]
  [, ..bv_keep_fields]
)
bv_plot_data <- merge(
  x = adm_boundaries,
  y = bv_data_dt,
  by = adm_ids
)

bv_fig <- ggplot() + 
  geom_sf(data = bv_plot_data, aes(fill = idx), linewidth = 0.25, color = '#444444') + 
  scale_fill_manual(values = bv_color_vec) +
  guides(fill = 'none') +
  theme_void()

full_fig <- bv_fig +
  (patchwork::plot_spacer() / color_scheme_fig / patchwork::plot_spacer()) +
  patchwork::plot_layout(nrow = 1, widths = c(6, 2))

pdf(file.path(viz_dir, 'bivariate_uncertainty_plot.pdf'), height = 8, width = 8)
plot(full_fig)
dev.off()

# Plot gridded population density across the country ------------------------------------>

pop_table <- raster_to_table(pop_raster, type = 'Population')
# Get reasonably-spaced breaks on a log scale
max_val <- max(pop_table$value, na.rm = T)
pow_10 <- floor(log(max_val) / log(10))
max_val_mult <- round(max_val / 10^pow_10) 
pop_breaks <- c(0, max_val_mult * 10^unique(seq(0, pow_10, length.out = 3)))

pop_fig <- ggplot() +
  geom_raster(data = pop_table, aes(x=x, y=y, fill=value)) +
  geom_sf(data = adm_boundaries, fill = NA, linewidth = 0.25, color = "black") + 
  scale_fill_gradientn(
    colors = viridisLite::viridis(10)[2:10], labels = scales::comma, breaks = pop_breaks,
    limits = range(pop_breaks), oob = scales::squish, trans = 'log1p'
  ) + 
  labs(title = "Population density", x = '', y = '', fill = "Population\n(log scaled)") +
  map_theme

pdf(file.path(viz_dir, 'population_density.pdf'), height = 8, width = 8)
plot(pop_fig)
dev.off()
