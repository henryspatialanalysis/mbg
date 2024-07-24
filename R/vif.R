#' Run variance inflation factor (VIF) selection on input covariates
#'
#' @param dataset data.frame-like object with named columns containing all covariates to
#'   consider in the VIF analysis.
#' @param vif_cutoff (`numeric(1)`) Cutoff for maximum variance inflation factor in
#'  `dataset`
#'
#' @return data.table listing each variable, VIF in most recent round, and whether the
#'   indicator should be included or not
#'
#' @import data.table glue
#' @importFrom stats lm
#' @export
vif_covariate_select <- function(dataset, vif_cutoff = 5){
  # Set dummy variables to avoid R CMD CHECK issues
  keep <- vif <- NULL
  
  # Check inputs
  if(!is.data.frame(dataset)) stop("VIF input dataset must be data.frame-like")
  if(max(table(colnames(dataset))) > 1) stop("VIF input dataset must have unique column names")
  if(!is.numeric(vif_cutoff)) stop("VIF cutoff must be numeric")
  if(is.na(vif_cutoff) | (vif_cutoff <= 0)) stop("VIF cutoff must be greater than zero")

  # Helper function to convert from R-squared to VIF
  r2_to_vif <- function(r2) 1 / (1 - r2)

  # Iteratively run VIF until we get to 2 columns or pass below the cutoff
  remaining <- colnames(dataset)
  if(length(remaining) <= 1L){
    return(data.table::data.table(covariate = remaining, vif = NA_real_, keep = TRUE))
  }

  results_list <- list()
  run_again <- TRUE
  while(run_again){
    vifs_this_round <- sapply(remaining, function(check_col){
      other_cols <- setdiff(remaining, check_col)
      model <- stats::lm(
        glue::glue("{check_col} ~ {paste(other_cols, collapse = ' + ')}"),
        data = dataset
      )
      return(r2_to_vif(summary(model)$r.squared))
    })
    max_vif <- max(vifs_this_round)
    if((max_vif > vif_cutoff) & length(remaining) > 3){
      max_vif_name <- remaining[which.max(vifs_this_round)]
      message(glue::glue("VIF for {max_vif_name} = {round(max_vif, 4)} > {vif_cutoff}"))
      results_list <- c(
        results_list,
        list(data.table::data.table(covariate = max_vif_name, vif = max_vif))
      )
      remaining <- setdiff(remaining, max_vif_name)
    } else {
      results_list <- c(
        results_list,
        list(data.table::data.table(covariate = remaining, vif = vifs_this_round))
      )
      run_again <- FALSE
    }
  }
  # Combine results into a single table
  vif_results_table <- data.table::rbindlist(results_list)[, keep := vif < vif_cutoff ]
  return(vif_results_table)
}
