## #######################################################################################
##
## Validate MBG results by comparing to statcompiler data at the adm1 level
##
## AUTHOR: Nathaniel Henry, nat@henryspatialanalysis.com
## CREATED: 26 September 2023
## PURPOSE: Check outputs of mbg_package_workflow.R at the admin1 level
##
## NOTE: Relies on the `rdhs` and `fedmatch` packages. To install, run these lines:
## install.packages('fedmatch')
## devtools::install_github("ropensci/rdhs")
##
## #######################################################################################

## 00) SETTINGS ------------------------------------------------------------------------->

REPOS_PATH_DEFAULT <- '~/repos'

# The script takes one required argument:
#  - `config_path`: the path to the `config.yaml` configuration file that was saved in
#    mbg_package_workflow.R
#  - `statcompiler_id`: the indicator ID in StatCompiler
#  - `statcompiler_denom`: Does the StatCompiler Indicator need to be divided by a
#    denominator to match the MBG data (e.g., set as 100 to convert from percentage to a
#    true fraction). If NULL, does not run any conversion
#
# The script also takes one optional argument, which has a default:
#  - `repos_path`: the path to the directory that contains the `versioning` custom package

if(interactive()){
  config_path <- '~/temp_data/geostats/mbg_results/20230926/config.yaml'
  statcompiler_id <- 'CN_NUTS_C_WH2'
  statcompiler_denom <- 100
  repos_path <- REPOS_PATH_DEFAULT
} else {
  library(argparse)
  parser <- argparse::ArgumentParser()
  parser$add_argument("--config_path", type = 'character')
  parser$add_argument("--statcompiler_id", type = 'character')
  parser$add_argument("--statcompiler_denom", type = 'character', default = NULL)
  parser$add_argument("--repos_path", type = 'character', default = REPOS_PATH_DEFAULT)
  # Parse command line arguments
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  config_path <- args$config_path
  statcompiler_id <- args$statcompiler_id
  statcompiler_denom <- args$statcompiler_denom
  repos_path <- args$repos_path
}


## 01) SETUP ---------------------------------------------------------------------------->

# Load required packages
load_packages <- c('data.table', 'fedmatch', 'ggplot2', 'glue', 'rdhs', 'scales')
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
dir.create(viz_dir, showWarnings = FALSE)

# Load MBG results for visualization
adm1_mbg_summaries <- config$read("results", "adm1_summary_table")


## 02) Prepare StatCompiler data -------------------------------------------------------->

# Get the country code associated with this country
country_lookup <- dhs_countries(returnFields = c("CountryName", "DHS_CountryCode")) |>
  data.table::as.data.table()
country_id <- country_lookup[CountryName == config$get("country"), DHS_CountryCode]
if(length(country_id) != 1) stop("Did not return a unique country code - check name")  

# Pull data from DHS StatCompiler API
sc_summaries <- rdhs::dhs_data(
  countryIds = country_id,
  indicatorIds = statcompiler_id,
  surveyYearStart = config$get('year'),
  breakdown = 'subnational'
) |> as.data.table()
sc_summaries[, adm_name := gsub('\\.', '', CharacteristicLabel)]
sc_summaries[, dhs_adm_level := (nchar(CharacteristicLabel) - nchar(adm_name))/2 + 1]
if(!is.null(statcompiler_denom)) sc_summaries[, Value := Value / statcompiler_denom ]


## 03) Merge StatCompiler and MBG results using fuzzy string matching ------------------->

fuzzy_merge_results <- fedmatch::merge_plus(
  data1 = adm1_mbg_summaries[, .(ADM1_NAME, ADM1_CODE, mean, lower, upper)],
  data2 = sc_summaries[, .(adm_name, dhs_adm_level, sc_value = Value, sc_id = .I)],
  by.x = 'ADM1_NAME',
  by.y = 'adm_name',
  match_type = 'fuzzy',
  unique_key_1 = 'ADM1_CODE',
  unique_key_2 = 'sc_id'
)
fuzzy_matches <- fuzzy_merge_results$matches
fuzzy_matches[, in_ui := 'No' ][(sc_value >= lower) & (sc_value <= upper), in_ui := 'Yes']
n_total <- nrow(adm1_mbg_summaries)
n_matched <- nrow(fuzzy_matches)
n_in_ui <- sum(fuzzy_matches$in_ui == "Yes")

summary_text <- glue::glue(
  "Of {n_total} admin1 units in the MBG shapefile, {n_matched} ",
  "({scales::percent(n_matched/n_total)}) were matched to DHS StatCompiler estimates.\n",
  "Of {n_matched} matched units, {n_in_ui} ({scales::percent(n_in_ui / n_matched)}) ",
  "StatCompiler estimates fell within the MBG 95% Uncertainty Interval."
)
message(summary_text)


## 04) Plot admin1 comparison ----------------------------------------------------------->

ui_colors <- c(Yes = 'black', No = 'red')

fig <- ggplot(data = fuzzy_matches, aes(x = sc_value, y = mean, color = in_ui)) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey40', linetype = 3, linewidth = 1) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_label(
    data = function(x) x[in_ui == "No", ],
    aes(label = adm_name),
    nudge_x = diff(range(fuzzy_matches$sc_value)) / 75,
    hjust = 0, label.size = NA, fill = alpha("white", 0.6), show.legend = FALSE
  ) +
  scale_color_manual(values = ui_colors) + 
  scale_x_continuous(labels = scales::percent) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "StatCompiler Estimates",
    y = "MBG Estimates (Aggregated)",
    color = "In MBG UI?",
    subtitle = summary_text
  ) +
  theme_bw()

pdf(file.path(viz_dir, "statcompiler_vetting.pdf"), height = 8, width = 8)
print(fig)
dev.off()




# See the config.yaml file for the list of all results files that can be loaded
modeling_level <- config$get('shapefile_settings', 'modeling_level')
adm_boundaries <- config$read('results', 'adm_boundaries')
summary_table <- config$read('results', glue::glue('adm{modeling_level}_summary_table'))
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
  geom_sf(data = adm_boundaries, fill = NA, linewidth = 0.05, color = '#444444') + 
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
  scale_fill_gradientn(
    colors = color_scheme, limits = value_limits, oob = scales::squish,
    labels = scales::percent
  ) +
  labs(
    title = glue::glue("Summary estimates: {indicator} in {country}, {year}"),
    x = '', y = '', fill = indicator
  ) +
  theme_minimal() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
pdf(file.path(viz_dir, 'summary_admin_estimates.pdf'), height = 6, width = 10)
print(adm_fig)
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
