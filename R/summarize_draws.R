#' Summarize draws
#'
#' @description Helper function to summarize a matrix or data.frame of predictive draws
#'
#' @param draws A `matrix`, `data.frame`, or [data.table::data.table] of predictive draws.
#' @param id_fields (default NULL) Only considered for data.frame-like `draws`. What
#'   identifier fields in the data should be kept in the summary table and not included
#'   among the draw fields?
#' @param draw_fields (default NULL) Only considered for data.frame-like `draws`. What
#'   fields represent actual draws, as opposed to identifier fields or other metadata like
#'   population? If `NULL`, the default, automatically determines the draw fields as all
#'   columns not included in the `id_fields`.
#' @param ui_width (`numeric`, default 0.95) Size of the uncertainty interval width when
#'   calculating the upper and lower summary rasters
#' @param na.rm (`logical`, default TRUE) Should NA values be removed when calculating
#'   summaries across draws?
#'
#' @return A [data.table::data.table] containing at least the following fields:
#'   - The `id_fields`, if passed
#'   - "mean": Mean across predictive draws
#'   - "lower": Lower bound of the (X%) uncertainty interval
#'   - "upper": Upper bound of the (X%) uncertainty interval
#'   - "ui_width": "upper" - "lower"
#'
#' @examples
#' # Summarize a draws matrix
#' draws_matrix <- matrix(rnorm(200), nrow = 10)
#' summary_table_a <- summarize_draws(draws_matrix)
#' head(summary_table_a)
#'
#' # Summarize a draws data.table with location IDs
#' draws_table <- matrix(c(1:10, rnorm(200)), nrow = 10) |>
#'   data.table::as.data.table() |>
#'   data.table::setnames(c('location_id', paste0('draw_', 1:20)))
#' summary_table_b <- summarize_draws(draws_table, id_fields = 'location_id')
#' head(summary_table_b)
#'
#' @concept prediction
#'
#' @import data.table
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowQuantiles
#' @export
summarize_draws <- function(
  draws, id_fields = NULL, draw_fields = NULL, ui_width = 0.95, na.rm = TRUE
){
  if(inherits(draws, 'data.frame')){
    draws <- data.table::as.data.table(draws)
    if(is.null(draw_fields)) draw_fields <- setdiff(colnames(draws), id_fields)
    if(!is.null(id_fields)){
      ids_table <- draws[, id_fields, with = F]
    } else {
      ids_table <- NULL
    }
    draws_mat <- as.matrix(draws[, draw_fields, with = F])
  } else {
    draws_mat <- draws
    ids_table <- NULL
  }
  # Summarize as a table
  alpha <- (1 - ui_width) / 2
  summary_table <- data.table::data.table(
    mean = Matrix::rowMeans(draws_mat, na.rm = na.rm),
    lower = matrixStats::rowQuantiles(draws_mat, probs = alpha, na.rm = na.rm),
    upper = matrixStats::rowQuantiles(draws_mat, probs = 1 - alpha, na.rm = na.rm)
  )
  summary_table$ui_width <- summary_table$upper - summary_table$lower
  if(!is.null(ids_table)) summary_table <- cbind(ids_table, summary_table)
  return(summary_table)
}
