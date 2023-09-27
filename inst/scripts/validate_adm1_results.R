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

# Save the raw StatCompiler data and merged estimates
versioning::autowrite(
  sc_summaries,
  file = file.path(viz_dir, 'statcompiler_estimates.csv')
)
versioning::autowrite(
  fuzzy_matches,
  file = file.path(viz_dir, 'statcompiler_estimates_matched.csv')
)


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
