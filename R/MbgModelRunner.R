#' @title MBG model runner class
#'
#' @description
#' R6 class to run a full MBG model and make predictions.
#'
#' @details To see examples of this object, run \code{vignette('mbg')}
#'
#' @concept model_runner
#'
#' @import data.table
#' @importFrom R6 R6Class
#' @importFrom matrixStats rowQuantiles
#' @export
MbgModelRunner <- R6::R6Class(
  "MbgModelRunner",
  public = list(

    ## Fields set during initialization ------------------------------------------------->

    #' @field input_data ([data.table::data.table])\cr
    #' Table containing at least the following fields:\cr
    #'   * x (`numeric`) location longitude in decimal degrees\cr
    #'   * y (`numeric`) location latitude in decimal degrees\cr
    #'   * indicator (`integer`) The number of events in the population\cr
    #'   * samplesize (`integer`) The total population, denominator for `indicator`
    input_data = NULL,

    #' @field id_raster ([terra::SpatRaster])\cr
    #' raster showing the total area that will be predicted using this model.
    id_raster = NULL,

    #' @field covariate_rasters (`list()`)\cr
    #' A list containing all predictor covariates. Each covariate is a
    #' [terra::SpatRaster] object with the same extent and dimensions as
    #' `id_raster`.
    covariate_rasters = NULL,

    #' @field aggregation_table ([data.table::data.table])\cr
    #' A table created by [build_aggregation_table], used to link each grid cell to
    #' higher-level administrative units.
    aggregation_table = NULL,

    #' @field aggregation_levels (`list()`)\cr
    #' A named list: for each named item, the name is the label for that aggregation
    #' level, and the value is a character vector of all fields in the original polygons
    #' to be used for aggregation at that level.
    aggregation_levels = NULL,

    #' @field population_raster ([terra::SpatRaster])\cr
    #' A raster giving population for each grid cell, to be used for population-weighted
    #' aggregation from grid cells to polygon boundaries. Should have the same dimensions
    #' as `id_raster`. If no population raster is passed and the results are aggregated,
    #' aggregation will be by simple mean rather than population-weighted mean
    population_raster = NULL,

    #' @field admin_bounds ([sf::sf])\cr
    #'  Polygons showing the boundaries of administrative divisions within the study
    #'   region. Only required if `use_admin_effect` OR `stacking_use_admin_bounds` is
    #'   `TRUE`.
    admin_bounds = NULL,

    #' @field admin_bounds_id (`character`)\cr
    #'  Field containing unique identifiers for `admin_bounds`, if passed.
    admin_bounds_id = NULL,

    #' @field use_covariates (`logical(1)`)\cr
    #' Should covariate effects be included in the predictive model?
    use_covariates = NULL,

    #' @field use_gp (`logical(1)`)\cr
    #' Should a smoothed spatial surface be included in the predictive model?
    use_gp = NULL,

    #' @field use_admin_effect (`logical(1)`)\cr
    #' Should IID administrative-level effects be included in the predictive model?
    use_admin_effect = NULL,

    #' @field use_nugget (`logical(1)`)\cr
    #' Should an IID effect by pixel be included in the predictive model?
    use_nugget = NULL,

    #' @field use_stacking (`logical(1)`)\cr
    #' Should machine learning submodels be trained to relate the covariate rasters with
    #'   the outcome data? Only run if `use_covariates` is `TRUE`.
    use_stacking = NULL,

    #' @field stacking_model_settings (`list()`)\cr
    #' A named list of submodels to be run. For more information about this term, see
    #'   [run_regression_submodels]. Only considered if `use_stacking` is TRUE.
    stacking_model_settings = NULL,

    #' @field stacking_cv_settings (`list()`)\cr
    #' How should the stacking submodels be cross-validated? For more information about
    #'   this term, see [run_regression_submodels]. Only considered if `use_stacking` is
    #'  `TRUE`.
    stacking_cv_settings = NULL,

    #' @field stacking_use_admin_bounds (`logical(1)`)\cr
    #' Should admin boundaries be included as features in the stacking submodels? For more
    #'   information about this term, see [run_regression_submodels]. Only considered if
    #'   `use_stacking` is TRUE.
    stacking_use_admin_bounds = NULL,

    #' @field stacking_prediction_range (`logical(1)`)\cr
    #' Range of possible predictions for the stacking submodels. For more information
    #'   about this term, see [run_regression_submodels]. Only considered if
    #'   `use_stacking` is TRUE.
    stacking_prediction_range = NULL,

    #' @field mesh_max_edge (`numeric(2)` or NULL)\cr
    #' Maximum size of the INLA SPDE mesh inside (1) and outside (2) of the modeled
    #'   region. Only considered if `use_gp` is TRUE.
    mesh_max_edge = NULL,

    #' @field mesh_cutoff (`numeric(1)`)\cr
    #' Minimum size of the INLA mesh, usually reached in data-dense areas. Only considered
    #'   if `use_gp` is TRUE.
    mesh_cutoff = NULL,

    #' @field spde_integrate_to_zero (`boolean(1)`)\cr
    #' Should the 'volume' under the SPDE mesh integrate to zero? Only considered if
    #'   `use_gp` is TRUE.
    spde_integrate_to_zero = NULL,

    #' @field prior_spde_range (`list()`)\cr
    #' A named list specifying the penalized complexity prior for the SPDE range. The two
    #'   named items are "threshold", the test threshold (set as a proportion of the
    #'   overall mesh extent), and "prob_below", the prior probability that the value is
    #'   BELOW that range threshold. The function automatically converts "threshold" from
    #'   a proportion of the overall mesh extent into a distance. Only considered if
    #'   `use_gp` is TRUE.
    prior_spde_range = NULL,

    #' @field prior_spde_sigma (`list()`)\cr
    #' A named list specifying the penalized complexity prior for sigma (standard
    #'   deviation) of the SPDE object. The two named items are "threshold", the test
    #'   threshold for the standard deviation, and "prob_above", the prior probability
    #'   that sigma will EXCEED that threshold. Only considered if `use_gp` is TRUE
    prior_spde_sigma = NULL,

    #' @field prior_nugget (`list()`)\cr
    #' A named list specifying the penalized complexity prior for the nugget term. The two
    #'   named items are "threshold", the test threshold for the nugget standard
    #'   deviation, and "prob_above", the prior probability that the standard deviation
    #'   will EXCEED that threshold. Only considered if `use_nugget` is TRUE.
    prior_nugget = NULL,

    #' @field prior_admin_effect (`list()`)\cr
    #' A named list specifying the penalized complexity prior for the admin-level IID
    #'   term. The two named items are "threshold", the test threshold for the standard
    #'   deviation of admin-level effects, and "prob_above", the prior probability that
    #'   the standard deviation will EXCEED that threshold. Only considered if
    #'   `use_admin_effect` is TRUE.
    prior_admin_effect = NULL,

    #' @field prior_covariate_effect (`list()`)\cr
    #' A named list specifying the penalized complexity prior for all covariate effects
    #'   except for the intercept, if an intercept is included. The two named items are
    #'   "threshold", the test threshold for the size of each fixed effect, and
    #'   "prob_above", the prior probability that the beta for each covariate will EXCEED
    #'   that threshold. Only considered if `use_covariates` is TRUE and `use_stacking` is
    #'   FALSE.
    prior_covariate_effect = NULL,

    #' @field inla_link (`character(1)`)\cr
    #' Link function for fitting the INLA model, typically related to the GLM `family`.
    inla_link = NULL,

    #' @field inverse_link (`character(1)`)\cr
    #' Inverse function of `inla_link`.
    inverse_link = NULL,

    #' @field inla_family (character)\cr
    #' GLM family to use. For more information, see [stats::family].
    inla_family = NULL,

    #' @field nugget_in_predict (`logical(1)`)\cr
    #' If the nugget is used in model fitting, should it also be included as an IID effect
    #'   by pixel in the model prediction step?
    nugget_in_predict = NULL,

    #' @field verbose
    #' Should model progress be timed?
    verbose = TRUE,

    ## Fields set during model runs ----------------------------------------------------->

    #' @field model_covariates (`list()`)\cr
    #' A list of covariates to be included in the INLA model. Either equal to
    #'   `covariate_rasters`, or ML model predictions for stacked generalization.
    model_covariates = NULL,

    #' @field inla_inputs_list (`list()`)\cr
    #' List of model inputs yielded by [prepare_inla_data_stack]
    inla_inputs_list = NULL,

    #' @field inla_fitted_model (`list()`)\cr
    #' List of model outputs yielded by [fit_inla_model]
    inla_fitted_model = NULL,

    #' @field grid_cell_predictions
    #' List of predictive surfaces yielded by [generate_cell_draws_and_summarize]
    grid_cell_predictions = NULL,

    #' @field aggregated_predictions
    #' List of predictions by administrative unit. Only created if `aggregation_table` and
    #' `aggregation_levels` are both defined.
    aggregated_predictions = NULL,

    ## Public methods ------------------------------------------------------------------->

    #' @description
    #' Create a new MbgModelRunner object
    #'
    #' @param input_data ([data.table::data.table]) Table containing at least
    #'   the following fields:\cr
    #'   * x (`numeric`) location x position, in the same projection as the `id_raster`\cr
    #'   * y (`numeric`) location y position, in the same projection as the `id_raster`\cr
    #'   * indicator (`integer`) The number of events in the population\cr
    #'   * samplesize (`integer`) The total population, denominator for `indicator`\cr
    #' @param id_raster ([terra::SpatRaster]) raster showing the total area that
    #'   will be predicted using this model
    #' @param covariate_rasters (`list()`, default NULL) A list containing all predictor
    #'   covariates. Each covariate is a [terra::SpatRaster] object with the same extent
    #'   and dimensions as `id_raster`.
    #' @param aggregation_table ([data.table::data.table]) A table created by
    #'   [build_aggregation_table], linking each grid cell to one or more polygons
    #' @param aggregation_levels (`list()`) A named list: for each named item, the name is
    #' the label for that aggregation level, and the value is a character vector of all
    #' fields in the original polygons to be used for aggregation at that level.
    #' @param population_raster ([terra::SpatRaster]) A raster giving population for each
    #' grid cell, to be used for population-weighted aggregation from grid cells to
    #' polygon boundaries. Should have the same dimensions as `id_raster`. If no
    #' population raster is passed and the results are aggregated, aggregation will be by
    #' simple mean rather than population-weighted mean
    #' @param admin_bounds ([sf::sf], default `NULL`) Polygons showing the boundaries
    #'   of administrative divisions within the study region. Only required if
    #'   `use_admin_effect` OR `stacking_use_admin_bounds` is `TRUE`.
    #' @param admin_bounds_id (`character`, default `NULL`) Field containing unique
    #'   identifiers for `admin_bounds`, if passed.
    #' @param use_covariates (`logical(1)`, default TRUE) Should covariate effects be
    #'   included in the predictive model?
    #' @param use_gp (`logical(1)`, default TRUE) Should a smoothed spatial surface be
    #'   included in the predictive model?
    #' @param use_admin_effect (`logical(1)` default FALSE) Should IID
    #'   administrative-level effects be included in the predictive model?
    #' @param use_nugget (`logical(1)`, default TRUE) Should an IID effect by pixel be
    #'   included in the predictive model?
    #' @param use_stacking (`logical(1)`, default FALSE) Should machine learning submodels
    #'   be trained to relate the covariate rasters with the outcome data? Only run if
    #'   `use_covariates` is `TRUE`.
    #' @param stacking_model_settings (`list()`) A named list of submodels to be run. For
    #'   more information about this term, see [run_regression_submodels]. Only considered
    #'   if `use_stacking` is TRUE.
    #' @param stacking_cv_settings (`list()`) How should the stacking submodels be
    #'   cross-validated? For more information about this term, see
    #'  [run_regression_submodels]. Only considered if `use_stacking` is `TRUE`.
    #' @param stacking_use_admin_bounds (`logical(1)`, default FALSE) Should admin
    #'   boundaries be included as features in the stacking submodels? For more
    #'   information about this term, see [run_regression_submodels]. Only considered if
    #'   `use_stacking` is TRUE.
    #' @param stacking_prediction_range (`numeric(2)`, default NULL) Range of possible
    #'   predictions for the stacking submodels. For more information about this term, see
    #'   [run_regression_submodels]. Only considered if `use_stacking` is TRUE.
    #' @param mesh_max_edge (`numeric(2)`, default c(0.2, 5)) Maximum size of the INLA
    #'   SPDE mesh inside (1) and outside (2) of the modeled region. Only considered if
    #'   `use_gp` is TRUE.
    #' @param mesh_cutoff (`numeric(1)`, default 0.04) Minimum size of the INLA mesh,
    #'   usually reached in data-dense areas. Only considered if `use_gp` is TRUE.
    #' @param spde_integrate_to_zero (`boolean(1)`, default FALSE) Should the 'volume'
    #'   under the SPDE mesh integrate to zero? Only considered if `use_gp` is TRUE.
    #' @param prior_spde_range (`list()`) A named list specifying the penalized complexity
    #'   prior for the SPDE range. The two named items are "threshold", the test threshold
    #'   (set as a proportion of the overall mesh extent), and "prob_below", the prior
    #'   probability that the value is BELOW that range threshold. The function
    #'   automatically converts "threshold" from a proportion of the overall mesh extent
    #'   into a distance. Only considered if `use_gp` is TRUE.
    #' @param prior_spde_sigma (`list()`) A named list specifying the penalized complexity
    #'   prior for sigma (standard deviation) of the SPDE object. The two named items are
    #'   "threshold", the test threshold for the standard deviation, and "prob_above",
    #'   the prior probability that sigma will EXCEED that threshold. Only considered if
    #'   `use_gp` is TRUE
    #' @param prior_nugget (`list()`) A named list specifying the penalized complexity
    #'   prior for the nugget term. The two named items are "threshold", the test
    #'   threshold for the nugget standard deviation, and "prob_above", the prior
    #'   probability that the standard deviation will EXCEED that threshold. Only
    #'   considered if `use_nugget` is TRUE.
    #' @param prior_admin_effect (`list()`) A named list specifying the penalized
    #'   complexity prior for the admin-level IID term. The two named items are
    #'   "threshold", the test threshold for the standard deviation of admin-level
    #'   effects, and "prob_above", the prior probability that the standard deviation will
    #'   EXCEED that threshold. Only considered if `use_admin_effect` is TRUE.
    #' @param prior_covariate_effect (`list()`) A named list specifying the penalized
    #'   complexity prior for all covariate effects except for the intercept, if an
    #'   intercept is included. The two named items are "threshold", the test threshold
    #'   for the size of each fixed effect, and "prob_above", the prior probability that
    #'   the beta for each covariate will exceed that threshold. Only considered if
    #'   `use_covariates` is TRUE and `use_stacking` is FALSE.
    #' @param inla_link (`character(1)`, default 'logit') Link function for fitting the
    #'   INLA model, typically related to the GLM `family`.
    #' @param inverse_link (`character(1)`, default 'plogis') Inverse function of
    #'   `inla_link`.
    #' @param inla_family (`character(1)`, default 'binomial') GLM family to use. For more
    #'   information, see [stats::family()].
    #' @param nugget_in_predict (`logical(1)`, default TRUE) If the nugget is used in
    #'   model fitting, should it also be included as an IID effect by pixel in the
    #'   model prediction step?
    #' @param verbose (`logical(1)`, default TRUE) Should model progress be timed?
    #'
    initialize = function(
      input_data,
      id_raster,
      covariate_rasters = NULL,
      aggregation_table = NULL,
      aggregation_levels = NULL,
      population_raster = NULL,
      admin_bounds = NULL,
      admin_bounds_id = NULL,
      use_covariates = TRUE,
      use_gp = TRUE,
      use_admin_effect = FALSE,
      use_nugget = TRUE,
      use_stacking = FALSE,
      stacking_cv_settings = list(
        method = 'repeatedcv',
        number = 5,
        repeats = 5
      ),
      stacking_model_settings = list(
        gbm = NULL,
        treebag = NULL,
        rf = NULL
      ),
      stacking_use_admin_bounds = FALSE,
      stacking_prediction_range = NULL,
      mesh_max_edge = c(0.2, 5.0),
      mesh_cutoff = c(0.04),
      spde_integrate_to_zero = FALSE,
      prior_spde_range = list(threshold = 0.1, prob_below = 0.05),
      prior_spde_sigma = list(threshold = 3, prob_above = 0.05),
      prior_nugget = list(threshold = 3, prob_above = 0.05),
      prior_admin_effect = list(threshold = 3, prob_above = 0.05),
      prior_covariate_effect = list(threshold = 3, prob_above = 0.05),
      inla_link = 'logit',
      inverse_link = 'plogis',
      inla_family = 'binomial',
      nugget_in_predict = TRUE,
      verbose = TRUE
    ){
      # Set values
      self$input_data <- input_data
      self$id_raster <- id_raster
      self$covariate_rasters <- covariate_rasters
      self$aggregation_table <- aggregation_table
      self$aggregation_levels <- aggregation_levels
      self$population_raster <- population_raster
      self$admin_bounds <- admin_bounds
      self$admin_bounds_id <- admin_bounds_id
      self$use_covariates <- use_covariates
      self$use_gp <- use_gp
      self$use_admin_effect <- use_admin_effect
      self$use_nugget <- use_nugget
      self$use_stacking <- use_stacking
      self$stacking_model_settings <- stacking_model_settings
      self$stacking_cv_settings <- stacking_cv_settings
      self$stacking_use_admin_bounds <- stacking_use_admin_bounds
      self$stacking_prediction_range <- stacking_prediction_range
      self$mesh_max_edge <- mesh_max_edge
      self$mesh_cutoff <- mesh_cutoff
      self$spde_integrate_to_zero <- spde_integrate_to_zero
      self$prior_spde_range <- prior_spde_range
      self$prior_spde_sigma <- prior_spde_sigma
      self$prior_nugget <- prior_nugget
      self$prior_admin_effect <- prior_admin_effect
      self$prior_covariate_effect <- prior_covariate_effect
      self$inla_link <- inla_link
      self$inverse_link <- inverse_link
      self$inla_family <- inla_family
      self$nugget_in_predict <- nugget_in_predict
      self$verbose <- verbose
      return(self)
    },

    #' @description Prepare covariates for MBG model fitting
    prepare_covariates = function(){
      # Optionally run stacking
      if(self$use_covariates & self$use_stacking){
        if(is.null(self$stacking_prediction_range)) stop("Set stacking_prediction_range")
        stackers_list <- run_regression_submodels(
          input_data = copy(self$input_data),
          id_raster = self$id_raster,
          covariates = self$covariate_rasters,
          cv_settings = self$stacking_cv_settings,
          model_settings = self$stacking_model_settings,
          family = self$inla_family,
          use_admin_bounds = self$stacking_use_admin_bounds,
          admin_bounds = self$admin_bounds,
          admin_bounds_id = self$admin_bounds_id,
          prediction_range = self$stacking_prediction_range,
          verbose = self$verbose
        )
        self$model_covariates <- stackers_list$predictions
      } else if(self$use_covariates){
        self$model_covariates <- self$covariate_rasters
      } else {
        self$model_covariates <- NULL
      }
      invisible(NULL)
    },

    #' @description Fit MBG model
    fit_mbg_model = function(){
      # Prepare data for MBG model run
      self$inla_inputs_list <- prepare_inla_data_stack(
        input_data = self$input_data,
        id_raster = self$id_raster,
        covariates = self$model_covariates,
        use_covariates = self$use_covariates,
        covariates_sum_to_one = self$use_stacking,
        family = self$inla_family,
        use_spde = self$use_gp,
        spde_range_pc_prior = self$prior_spde_range,
        spde_sigma_pc_prior = self$prior_spde_sigma,
        spde_integrate_to_zero = self$spde_integrate_to_zero,
        mesh_max_edge = self$mesh_max_edge,
        mesh_cutoff = self$mesh_cutoff,
        use_nugget = self$use_nugget,
        nugget_pc_prior = self$prior_nugget,
        use_admin_effect = self$use_admin_effect,
        admin_boundaries = self$admin_bounds,
        admin_pc_prior = self$prior_admin_effect
      )
      # Fit MBG model
      self$inla_fitted_model <- fit_inla_model(
        formula = self$inla_inputs_list$formula_string,
        data_stack = self$inla_inputs_list$inla_data_stack,
        spde = self$inla_inputs_list$spde,
        samplesize_vec = self$input_data$samplesize,
        precision_vec = (self$input_data$sd)**(-2),
        family = self$inla_family,
        link = self$inla_link,
        fixed_effects_pc_prior = self$prior_covariate_effect,
        verbose = self$verbose
      )
      invisible(NULL)
    },


    #' @description Generate predictions by grid cell
    #'
    #' @param n_samples (`integer(1)`, default 1000) Number of posterior predictive
    #'   samples to generate from the fitted model object.
    #' @param ui_width (`numeric(1)`, default 0.95) Uncertainty interval width. This
    #'   method will create summary rasters for quantiles ((1 - ui_width)/2) and
    #'   (1 - (1 - ui_width)/2).
    generate_predictions = function(n_samples = 1e3, ui_width = 0.95){
      if(is.null(self$inla_fitted_model)){
        stop("Must fit geostatistical model before generating predictions")
      }
      self$grid_cell_predictions <- generate_cell_draws_and_summarize(
        inla_model = self$inla_fitted_model,
        inla_mesh = self$inla_inputs_list$mesh,
        n_samples = n_samples,
        id_raster = self$id_raster,
        covariates = self$model_covariates,
        inverse_link_function = self$inverse_link,
        nugget_in_predict = self$nugget_in_predict,
        admin_boundaries = self$admin_bounds,
        ui_width = ui_width,
        verbose = self$verbose
      )
      invisible(self$grid_cell_predictions)
    },

    #' @description Aggregate grid cell predictions
    #'
    #' @param ui_width (`numeric(1)`, default 0.95) Uncertainty interval width. This
    #'   method will create summary "upper" and "lower" fields rasters for quantiles
    #'   ((1 - ui_width)/2) and (1 - (1 - ui_width)/2).
    #'
    #' @return List with the same names as `self$aggregation_levels`, aggregating by the
    #'   columns specified in `self$aggregation_levels`
    aggregate_predictions = function(ui_width = 0.95){
      # Validate inputs
      if(is.null(self$grid_cell_predictions)){
        stop("Must run generate_predictions() before aggregate_predictions()")
      }
      if(is.null(self$aggregation_table) || is.null(self$aggregation_levels)){
        stop(
          "Both aggregation_table and aggregation_levels must be set before running ",
          "aggregate_predictions()"
        )
      }
      # Aggregate by population-weighted mean if population raster is available, or simple
      #  mean otherwise
      if(is.null(self$population_raster)){
        message("Running unweighted aggregation from grid cells to polygons")
        agg_method <- 'mean'
      } else {
        message("Running population-weighted aggregation from grid cells to polygons")
        agg_method <- 'weighted.mean'
      }
      # Run aggregation at all specified levels
      self$aggregated_predictions <- lapply(self$aggregation_levels, function(cols){
        aggregated_draws <- aggregate_draws_to_polygons(
          draws_matrix = self$grid_cell_predictions$cell_draws,
          aggregation_table = self$aggregation_table,
          aggregation_cols = cols,
          method = agg_method,
          weighting_raster = self$population_raster
        )
        # Summarize using mean and boundaries of the uncertainty interval
        draw_fields <- setdiff(colnames(aggregated_draws), cols)
        agg_draws_matrix <- as.matrix(aggregated_draws[, draw_fields, with = FALSE])
        alpha <- (1 - ui_width) / 2
        aggregated_summary <- data.table::copy(aggregated_draws[, cols, with = FALSE])
        aggregated_summary$mean <- rowMeans(agg_draws_matrix, na.rm = TRUE)
        aggregated_summary$lower <- matrixStats::rowQuantiles(
          agg_draws_matrix,
          probs = alpha,
          na.rm = TRUE
        )
        aggregated_summary$upper <- matrixStats::rowQuantiles(
          agg_draws_matrix,
          probs = 1 - alpha,
          na.rm = TRUE
        )
        # Add grid-cell population aggregated by polygon
        admin_pop <- aggregate_raster_to_polygons(
          data_raster = self$population_raster,
          aggregation_table = self$aggregation_table,
          aggregation_cols = cols,
          method = 'sum',
          aggregated_field = 'population'
        )
        aggregated_summary[admin_pop, population := i.population, on = cols]
        return(list(draws = aggregated_draws, summary = aggregated_summary))
      })
      names(self$aggregated_predictions) <- names(self$aggregation_levels)
    },

    #' @description
    #' Run a full MBG pipeline, including stacking, MBG model fitting, and prediction
    #'
    #' @param n_samples (`integer(1)`, default 1000) Number of posterior predictive
    #'   samples to generate from the fitted model object.
    #' @param ui_width (`numeric(1)`, default 0.95)
    run_mbg_pipeline = function(n_samples = 1e3, ui_width = 0.95){
      self$prepare_covariates()
      self$fit_mbg_model()
      self$generate_predictions(n_samples = n_samples, ui_width = ui_width)
      if(!is.null(self$aggregation_table) && !is.null(self$aggregation_levels)){
        self$aggregate_predictions()
        invisible(self$aggregated_predictions)
      } else {
        invisible(self$grid_cell_predictions)
      }
    },

    #' @description
    #' Get predictive validity metrics for the fitted model
    #'
    #' @details Returns the point RMSE (compared against the mean estimates by pixel),
    #'   log-posterior density (compared against the predictive draws), and the
    #'   Watanabe-Aikake Information Criterion (WAIC, only returned for in-sample
    #'   predictive validity).
    #'
    #' @param in_sample (`logical(1)`, default TRUE) Compare model predictions to
    #'   the data used to generate the model? If FALSE, does not return the WAIC, which
    #'   is only useful for in-sample predictive validity.
    #' @param validation_data ([data.table::data.table], default NULL) Observed data to
    #'   compare against. Expected for out-of-sample model validation. Table containing at
    #'   least the following fields:\cr
    #'   * x (`numeric`) location x position, in the same projection as the `id_raster`\cr
    #'   * y (`numeric`) location y position, in the same projection as the `id_raster`\cr
    #'   * indicator (`integer`) The number of events in the population\cr
    #'   * samplesize (`integer`) The total population, denominator for `indicator`
    #' @param na.rm (`logical(1)`, default FALSE) Should NA values be dropped from the
    #'   RMSE and log predictive density calculations?
    #'
    #' @return [data.table::data.table] Containing the following fields:\cr
    #'   * 'rmse': Root mean squared error when compared against the mean estimates by
    #'     pixel. Lower RMSE is better.\cr
    #'   * 'lpd': Log posterior predictive density when compared against pixel-level
    #'     samples from the model. Higher LPD is better.\cr
    #'   * 'waic' (in-sample only): Watanable-Aikake information criterion estimated by
    #'     INLA. Lower WAIC is better.\cr
    #' For clarity, these fields will have the suffix "_is" for in-sample models, and
    #' "_oos" for out-of-sample models.
    get_predictive_validity = function(
      in_sample = TRUE,
      validation_data = NULL,
      na.rm = FALSE
    ){
      # Check inputs
      if(is.null(self$grid_cell_predictions)){
        stop("Must run generate_predictions() before measure_predictive_validity()")
      }
      if(in_sample && !is.null(validation_data)){
        warning(paste(
          "Measuring in-sample predictive validity, but an external data source was",
          "provided. Are you sure you are not measuring out-of-sample?"
        ))
      }
      if(!in_sample && is.null(validation_data)){
        warning(paste(
          "Measuring out-of-sample predictive validity, but no external data source was",
          "provided. Model input data will be used."
        ))
      }
      if(is.null(validation_data)) validation_data <- data.table::copy(self$input_data)
      validation_data[, data_rate := indicator / samplesize]
      pv_table <- data.table::data.table(
        rmse = rmse_raster_to_point(
          estimates = self$grid_cell_predictions$cell_pred_mean,
          validation_data = validation_data,
          outcome_field = 'data_rate',
          na.rm = na.rm
        ),
        lpd = log_posterior_density(
          draws = self$grid_cell_predictions$cell_draws,
          validation_data = validation_data,
          id_raster = self$id_raster,
          na.rm = na.rm
        )
      )
      if(in_sample){
        # WAIC is returned only for in-sample predictive validity
        pv_table$waic <- self$inla_fitted_model$waic$waic
        # Append "_is" to all metrics
        colnames(pv_table) <- paste0(colnames(pv_table), "_is")
      } else {
        # Append "_oos" to all metrics
        colnames(pv_table) <- paste0(colnames(pv_table), "_oos")
      }
      return(pv_table)
    }
  )
)
