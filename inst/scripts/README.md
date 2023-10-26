# Scripts for `mbg`

The following scripts can be used alongside the `mbg` package functions:

- `mbg_package_workflow.R`: Run a geostatistics modeling workflow. This script relies on the partner packages `versioning` and `pixel2poly`, as well as a configuration file such as `inst/extdata/config.yaml` and a covariate table such as `inst/extdata/example_covariates.csv`.
- `visualize_mbg_results.R`: Visualize the results from `mbg_package_workflow.R`
- `validate_adm1_results.R`: Compare the estimates from `mbg_package_workflow.R`, aggregated at the admin1 level, to subnational estimates from the DHS StatCompiler tool. This script requires information about the country and indicator in the StatCompiler format.
