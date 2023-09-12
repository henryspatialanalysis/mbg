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
  config_path <- '~/temp_data/geostats/mbg_results/20230911/config.yaml'
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
adm2_boundaries <- config$read('results', 'adm2_shp')
adm2_summary_table <- config$read('results', 'adm2_summary_table')
mean_raster <- config$read('results', 'cell_pred_mean')
lower_raster <- config$read('results', 'cell_pred_lower')
upper_raster <- config$read('results', 'cell_pred_upper')
mesh <- config$read('results', 'inla_data_stack')$mesh

# Set some global variables that will be used in labels several times
indicator <- config$get('indicator')
country <- config$get('country')
year <- config$get('year')

# Set the color scheme for mean/lower/upper values
low_is_better <- TRUE
color_scheme <- RColorBrewer::brewer.pal(n = 9, name = 'Spectral')
if(low_is_better) color_scheme <- rev(color_scheme)


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
  geom_sf(data = adm2_boundaries, fill = NA, linewidth = 0.05, color = '#444444') + 
  scale_fill_gradientn(
    colors = color_scheme, limits = value_limits, oob = scales::squish,
    labels = scales::percent
  ) +
  labs(
    title = glue::glue("Summary rasters: {indicator} in {country}, {year}"),
    x = '', y = '', fill = indicator
  ) +
  theme_minimal() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
pdf(file.path(viz_dir, 'summary_rasters.pdf'), height = 6, width = 10)
print(raster_fig)
dev.off()


## 03) Plot the summary administrative data --------------------------------------------->

adm2_data_for_plotting <- merge(
  x = adm2_boundaries,
  y = melt(adm2_summary_table, id.vars = config$get("admin_id_fields")),
  by = config$get("admin_id_fields")
) |> merge(
  y = data.frame(
    variable = c('mean','lower','upper'),
    type = c('1. Mean', '2. Lower', '3. Upper')
  ),
  by = 'variable'
)

adm2_fig <- ggplot() + 
  facet_wrap('type', nrow = 1) +
  geom_sf(
    data = adm2_data_for_plotting,
    aes(fill = value), linewidth = 0.25, color = '#444444'
  ) + 
  scale_fill_gradientn(
    colors = color_scheme, limits = value_limits, oob = scales::squish,
    labels = scales::percent
  ) +
  labs(
    title = glue::glue("Summary district estimates: {indicator} in {country}, {year}"),
    x = '', y = '', fill = indicator
  ) +
  theme_minimal() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
pdf(file.path(viz_dir, 'summary_admin2_estimates.pdf'), height = 6, width = 10)
print(adm2_fig)
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

bv_keep_fields <- c(config$get('admin_id_fields'), 'idx')
bv_data_dt <- (copy(adm2_summary_table)
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
  x = adm2_boundaries,
  y = bv_data_dt,
  by = config$get('admin_id_fields')
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
