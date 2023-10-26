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
load_packages <- c(
  'data.table', 'fedmatch', 'ggplot2', 'glue', 'rdhs', 'scales', 'patchwork', 'cowplot'
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

# Save the raw StatCompiler data and merged estimates
versioning::autowrite(
  sc_summaries,
  file = file.path(viz_dir, 'statcompiler_estimates.csv')
)
versioning::autowrite(
  fuzzy_matches,
  file = file.path(viz_dir, 'statcompiler_estimates_matched.csv')
)

# Run a simple OLS model to see the relationship between MBG mean predictions and
#  StatCompiler estimates
ols_model <- stats::lm(mean ~ sc_value, data = fuzzy_matches)
ols_r_squared <- summary(ols_model)$r.squared |> round(3)
ols_coefficients <- summary(ols_model)$coefficients[, 1] |> round(3) |> unname()


## 04) Plot admin1 comparison ----------------------------------------------------------->

ui_colors <- c(Yes = '#444444', No = 'red')

subtitle_text <- glue::glue(
  '{summary_text}\nOLS regression: (MBG mean est.) ~ {ols_coefficients[1]} + ',
  '{ols_coefficients[2]} * (StatCompiler est.), R\u00b2 = {ols_r_squared}'
)

fig <- ggplot(data = fuzzy_matches, aes(x = sc_value, y = mean, color = in_ui)) + 
  geom_smooth(method = 'lm', col = '#339966', fill = '#99ffcc', alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = 'grey40', linetype = 3, linewidth = 1) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_label(
    data = function(x) x[in_ui == "No", ],
    aes(label = adm_name),
    nudge_x = diff(range(fuzzy_matches$sc_value)) / 75,
    hjust = 0, label.size = NA, fill = alpha("white", 0.8), show.legend = FALSE
  ) +
  scale_color_manual(values = ui_colors) + 
  scale_x_continuous(labels = scales::percent, limits = c(min(fuzzy_matches$sc_value) - 0.02, max(fuzzy_matches$sc_value) + .06)) + 
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "StatCompiler Estimates",
    y = "MBG Estimates (Aggregated)",
    color = "In MBG UI?",
    subtitle = subtitle_text
  ) +
  theme_bw()

pdf(file.path(viz_dir, "statcompiler_admin1_scatter.pdf"), height = 8, width = 8)
print(fig)
dev.off()


## 05) Plot MBG vs. StatCompiler absolute and relative differences ---------------------->

fuzzy_matches[, abs_diff := mean - sc_value ][, rel_diff := mean / sc_value ]
n_bins <- 7

# Extract a legend for the UI fill colors
for_legend <- ggplot(data = fuzzy_matches) +
  geom_histogram(aes(abs_diff, fill = in_ui)) + 
  scale_fill_manual(values = ui_colors) +
  labs(fill = "In MBG UI?")
diff_plot_legend <- cowplot::get_legend(for_legend)

# Build the absolute differences plot
abs_diff_mean <- fuzzy_matches[, round(mean(abs_diff), 3) ]
abs_diff_median <- fuzzy_matches[, round(median(abs_diff), 3) ]
abs_diff_plot <- ggplot(data = fuzzy_matches) + 
  geom_histogram(aes(abs_diff, fill = in_ui), center = 0, color = 'black', bins = n_bins) + 
  geom_vline(xintercept = 0, color = 'grey40', linetype = 3, linewidth = 1) +
  scale_fill_manual(values = ui_colors, guide = "none") +
  labs(
    title = 'Absolute differences',
    subtitle = glue::glue("Mean: {abs_diff_mean}; Median: {abs_diff_median}"),
    x = 'MBG mean est. - StatCompiler est.',
    y = 'Admin1 count'
  ) +
  theme_bw() +
  theme(axis.title = element_text(face = 'italic'))

# Build the relative differences plot
rel_diff_mean <- scales::percent(fuzzy_matches[, mean(rel_diff) ], accuracy = 1)
rel_diff_median <- scales::percent(fuzzy_matches[, median(rel_diff) ], accuracy = 1)
rel_diff_plot <- ggplot(data = fuzzy_matches) + 
  geom_histogram(aes(rel_diff, fill = in_ui), center = 1, color = 'black', bins = n_bins) + 
  geom_vline(xintercept = 1, color = 'grey40', linetype = 3, linewidth = 1) +
  scale_fill_manual(values = ui_colors, guide = "none") +
  scale_x_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = n_bins)) +
  labs(
    title = 'Relative differences',
    subtitle = glue::glue("Mean: {rel_diff_mean}; Median: {rel_diff_median}"),
    x = 'MBG mean est. / StatCompiler est.',
    y = 'Admin1 count'
  ) +
  theme_bw() +
  theme(axis.title = element_text(face = 'italic'))

## Combine plots and save

full_diff_plot <- ((abs_diff_plot / rel_diff_plot) | diff_plot_legend ) +
  patchwork::plot_layout(ncol = 2, widths = c(5, 1)) +
  patchwork::plot_annotation(
    title = 'Admin1 differences between MBG and StatCompiler',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

pdf(file.path(viz_dir, 'statcompiler_differences_histogram.pdf'), height = 8, width = 8)
plot(full_diff_plot)
dev.off()
