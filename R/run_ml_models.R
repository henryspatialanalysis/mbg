## #######################################################################################
##
## RUN REGRESSION SUB-MODELS 
##
## PURPOSE: Run regression models using the `caret` package for use in stacking
##
## #######################################################################################

#' Run regression sub-models
#' 
#' @description Wrapper to run many regression sub-models using the caret package
#' 
#' @param input_data
#' @param id_raster
#' @param covariates
#' @param cv_settings Cross-validation settings, passed to [caret::trainControl]
#' @param model_settings
#' 
#' @return List with two items:
#'   - "results": A list containing summary objects for each regression model
#'   - "predictions": Model predictions covering the entire id_raster
#'
#' @importFrom caret trainControl train predict
#' @importFrom terra extract
#' @import data.table
#' @export 
run_regression_submodels <- function(
  input_data, id_raster, covariates, cv_settings, model_settings
){
  # Prepare training data and eventual prediction space
  id_raster_table <- data.table::as.data.table(id_raster, xy = TRUE) |> na.omit()
  cov_names <- names(covariates)
  for(cov_name in cov_names){
    input_data[[cov_name]] <- terra::extract(
      x = covariates[[cov_name]],
      y = as.matrix(input_data[, .(x, y)])
    )[, 1]
    id_raster_table[[cov_name]] <- terra::extract(
      x = covariates[[cov_name]],
      y = as.matrix(id_raster_table[, .(x, y)])
    )[, 1]    
  }

  # Subset only to data outcome (indicator / samplesize), covariates, and x/y
  cov_cols <- c(cov_names, 'x', 'y')
  input_data[, data_rate := indicator / samplesize ]
  training_data <- input_data[, c('data_rate', cov_cols), with = F ]
  prediction_grid <- id_raster_table[, ..cov_cols ]

  # Run all prediction models
  fitControl <- do.call(trainControl, args = cv_settings)

  # Prepare prediction table
  model_types <- names(model_settings)
  models_list <- lapply(model_types, function(model_name){
    caret::train(
      "data_rate ~ .",
      data = training_data,
      method = model_name,
      model_settings[[model_name]]
    )
  })
  pred_list <- lapply(model_types, function(model_name){
    caret::predict(models_list[[model_name]], newdata = prediction_grid)
  })
  return(list(results = models_list, predictions = pred_list))
}