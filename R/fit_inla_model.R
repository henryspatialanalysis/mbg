#' Fit INLA model
#' 
#' @description Fit an INLA model based on a constructed data stack and formula
#' 
#' @details Using [INLA::inla()] with reasonable defaults and settings tuned to predict
#'   across a grid. NOTE: The sample size field in the `data_stack` MUST be named
#'   "samplesize".
#' 
#' @param formula INLA formula to fit. Generated in [prepare_inla_data_stack()]
#' @param data_stack Stacked data, covariates, and spatial index. Generated in
#'   [prepare_inla_data_stack()].
#' @param family (character, default 'binomial') GLM family to use. For more information,
#'   see [stats::family()].
#' @param link (character, default 'logit') Link function to use, typically related to the
#'   GLM `family`.
#' 
#' @return A fitted INLA model object created by [INLA::inla()]
#' 
#' @importFrom INLA inla inla.stack.data inla.stack.A
#' @importFrom tictoc tic toc
#' @export
fit_inla_model <- function(formula, data_stack, family = 'binomial', link = 'logit'){
  tictoc::tic("MBG model fitting")
  inla_model <- INLA::inla(
    formula = formula,
    family = family,
    Ntrials = samplesize,
    control.family = list(link = link),
    data = INLA::inla.stack.data(data_stack),
    control.predictor = list(
      compute = FALSE,
      link = 1,
      A = INLA::inla.stack.A(data_stack)
    ),
    control.compute = list(
      config = TRUE
    )
  )
  tictoc::toc()
  return(inla_model)
}
