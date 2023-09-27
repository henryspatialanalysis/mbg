#' Summarize draws
#' 
#' @description Helper function to summarize a matrix or data.frame of predictive draws
#' 
#' @param draws A matrix, data.frame, or data.table of predictive draws.
#' @param id_fields (default NULL) Only considered for data.frame-like `draws`. What
#'   identifier fields in the data should be kept in the summary table and not included
#'   among the draw fields?
#' @param draw_fields (default NULL) Only considered for data.frame-like `draws`. What 
#'   fields represent actual draws, as opposed to identifier fields or other metadata like
#'   population? If `NULL`, the default, automatically determines the draw fields as all
#'   columns not included in the `id_fields`.
#' @param ui_width (numeric, default 0.95) Size of the uncertainty interval width when
#'   calculating the upper and lower summary rasters
#' @param na.rm (logical, default TRUE) Should NA values be removed when calculating
#'   summaries across draws?
#' 
#' @return A data.table containing at least the following fields:
#'   - The `id_fields`, if passed
#'   - "mean": Mean across predictive draws
#'   - "lower": Lower bound of the (X%) uncertainty interval
#'   - "upper": Upper bound of the (X%) uncertainty interval
#' 
#' @import data.table
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowQuantiles
#' @export
summarize_draws <- function(
  draws, id_fields = NULL, draw_fields = NULL, ui_width = 0.95, na.rm = TRUE
){
  if(inherits(draws, 'data.frame')){
    draws <- as.data.table(draws)
    if(is.null(draw_fields)) draw_fields <- setdiff(colnames(draws), id_fields)
    if(!is.null(id_fields)){
      ids_table <- draws[, ..id_fields]
    } else {
      ids_table <- NULL
    }
    draws_mat <- as.matrix(draws[, ..draw_fields])
  } else {
    draws_mat <- draws
    ids_table <- NULL
  }
  # Summarize as a table
  summary_table <- data.table::data.table(
    mean = Matrix::rowMeans(draws_mat, na.rm = na.rm),
    lower = matrixStats::rowQuantiles(draws_mat, probs = (1 - ui_width)/2, na.rm = na.rm),
    upper = matrixStats::rowQuantiles(draws_mat, probs = 1 - (1 - ui_width)/2, na.rm = na.rm)
  )
  if(!is.null(ids_table)) summary_table <- cbind(ids_table, summary_table)
  return(summary_table)
}
