#' @title MBG model runner class
#'
#' @description
#' R6 class to run a full MBG model and make predictions.
#'
#' @importFrom R6 R6Class
#' @export
MbgModelRunner <- R6::R6Class(
  "MbgModelRunner",
  public = list(

    ## Fields set during initialization ------------------------------------------------->

    #' @field input_data ([data.table][data.table::data.table])\cr
    #' Table containing at least the following fields:\cr
    #'   * x (`numeric`) location longitude in decimal degrees
    #'   * y (`numeric`) location latitude in decimal degrees
    #'   * indicator (`integer`) The number of events in the population
    #'   * samplesize (`integer`) The total population, denominator for `indicator`
    input_data = NULL,

    #' @field id_raster ([terra][terra::rast])\cr
    #' raster showing the total area that will be predicted using this model.
    id_raster = NULL,

    #' @field covariate_rasters (`list()`)\cr
    #' A list containing all predictor covariates. Each covariate is a
    #' [terra][terra::rast] object with the same extent and dimensions as
    #' `id_raster`.
    covariate_rasters = NULL,

    #' @field admin_bounds ([sf][sf::sf])\cr
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
    #' GLM family to use. For more information, see [stats::family()].
    inla_family = NULL,

    #' @field nugget_in_predict (`logical(1)`)\cr
    #' If the nugget is used in model fitting, should it also be included as an IID effect
    #'   by pixel in the model prediction step?
    nugget_in_predict = NULL,


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


    ## Public methods ------------------------------------------------------------------->

    #' @description
    #' Create a new MbgModelRunner object
    #' 
    #' @param input_data ([data.table][data.table::data.table]) Table containing at least
    #'   the following fields:\cr
    #'   * x (`numeric`) location longitude in decimal degrees
    #'   * y (`numeric`) location latitude in decimal degrees
    #'   * indicator (`integer`) The number of events in the population
    #'   * samplesize (`integer`) The total population, denominator for `indicator`
    #' @param id_raster ([terra][terra::rast]) raster showing the total area that
    #'   will be predicted using this model
    #' @param covariate_rasters (`list()`) A list containing all predictor covariates.
    #'   Each covariate is a [terra][terra::rast] object with the same extent and
    #'   dimensions as `id_raster`.
    #' @param admin_bounds ([sf][sf::sf], default `NULL`) Polygons showing the boundaries
    #'   of administrative divisions within the study region. Only required if
    #'   `use_admin_effect` OR `stacking_use_admin_bounds` is `TRUE`.
    #' @param admin_bounds_id (`character`, default `NULL`) Field containing unique
    #'   identifiers for `admin_bounds`, if passed.
    #' @param use_covariates (`logical(1)`, default TRUE) Should covariate effects be
    #'   included in the predictive model?
    #' @param use_gp (`logical(1)`, default TRUE) Should a smoothed spatial surface be
    #'   included in the predictive model?
    #' @param use_admin_effect (`logical(1)` default TRUE) Should IID administrative-level
    #'   effects be included in the predictive model?
    #' @param use_nugget (`logical(1)`, default TRUE) Should an IID effect by pixel be
    #'   included in the predictive model?
    #' @param use_stacking (`logical(1)`, default TRUE) Should machine learning submodels
    #'   be trained to relate the covariate rasters with the outcome data? Only run if
    #'   `use_covariates` is `TRUE`.
    #' @param stacking_model_settings (`list()`) A named list of submodels to be run. For
    #'   more information about this term, see [run_regression_submodels]. Only considered
    #'   if `use_stacking` is TRUE.
    #' @param stacking_cv_settings (`list()`) How should the stacking submodels be
    #'   cross-validated? For more information about this term, see
    #'  [run_regression_submodels]. Only considered if `use_stacking` is `TRUE`.
    #' @param stacking_use_admin_bounds (`logical(1)`, default TRUE) Should admin
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
    #' @param spde_integrate_to_zero (`boolean(1)`, default TRUE) Should the 'volume'
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
    #'   the beta for each covariate will EXCEED that threshold. Only considered if
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
    #' 
    initialize = function(
      input_data,
      id_raster,
      covariate_rasters,
      admin_bounds = NULL,
      admin_bounds_id = NULL,
      use_covariates = TRUE,
      use_gp = TRUE,
      use_admin_effect = TRUE,
      use_nugget = TRUE,
      use_stacking = TRUE,
      stacking_cv_settings = list(
        method = 'repeatedcv',
        number = 5,
        repeats = 5
      ),
      stacking_model_settings = list(
        gam = NULL,
        gbm = NULL,
        treebag = NULL,
        rf = NULL
      ),
      stacking_use_admin_bounds = FALSE,
      stacking_prediction_range = NULL,
      mesh_max_edge = c(0.2, 5.0),
      mesh_cutoff = c(0.04),
      spde_integrate_to_zero = TRUE,
      prior_spde_range = list(threshold = 0.1, prob_below = 0.05),
      prior_spde_sigma = list(threshold = 3, prob_above = 0.05),
      prior_nugget = list(threshold = 3, prob_above = 0.05),
      prior_admin_effect = list(threshold = 3, prob_above = 0.05),
      prior_covariate_effect = list(threshold = 3, prob_above = 0.05),
      inla_link = 'logit',
      inverse_link = 'plogis',
      inla_family = 'binomial',
      nugget_in_predict = TRUE
    ){
      # Check internal consistency of inputs
      # TODO

      # Set values
      self$input_data <- input_data
      self$id_raster <- id_raster
      self$covariate_rasters <- covariate_rasters
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
      return(self)
    },

    #' @description Prepare covariates for MBG model fitting
    #' 
    #' @seealso [run_regression_submodels]
    prepare_covariates = function(){
      # Optionally run stacking
      if(self$use_covariates & self$use_stacking){
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
          prediction_range = self$stacking_prediction_range
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
    #' 
    #' @seealso [prepare_inla_data_stack], [fit_inla_model]
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
        fixed_effects_pc_prior = self$prior_covariate_effect
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
    #' 
    #' @seealso [generate_cell_draws_and_summarize]
    generate_predictions = function(n_samples = 1e3, ui_width = 0.95){
      self$grid_cell_predictions <- generate_cell_draws_and_summarize(
        inla_model = self$inla_fitted_model,
        inla_mesh = self$inla_inputs_list$mesh,
        n_samples = n_samples,
        id_raster = self$id_raster,
        covariates = self$model_covariates,
        inverse_link_function = self$inverse_link,
        nugget_in_predict = self$nugget_in_predict,
        admin_boundaries = self$admin_bounds,
        ui_width = ui_width
      )
      invisible(self$grid_cell_predictions)
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
      invisible(self$grid_cell_predictions)
    }
  )
)
