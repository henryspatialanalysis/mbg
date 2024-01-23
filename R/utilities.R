#' Make time stamp
#' 
#' @description Create a string time stamp based on current detailed date/time
#' 
#' @param suffix (`character(1)`, default NULL) suffix to append to the time stamp. Useful
#'   when running batches of related models 
#' @param milliseconds (`logical(1)`, default TRUE) Should milliseconds be appended to
#'   the timestamp? Useful when launching many models in quick succession.
#' 
#' @return A string formatted as 'YYYYMMDD_HH_MM_SS(_optional MS)(_optional suffix)'
#' 
#' @export
make_time_stamp <- function(suffix = NULL, milliseconds = T){
  if(milliseconds){
    time_stamp <- strftime(x = Sys.time(), format = '%Y%m%d_%H_%M_%OS3') |> 
      gsub(pattern = '\\.', replacement = '_')
  } else {
    time_stamp <- strftime(x = Sys.time(), format = '%Y%m%d_%H_%M_%S')
  }

  # Suffix must either be NULL or length 1
  if(length(suffix) > 1) stop("suffix should be NULL or a character vector of length 1.")
  if(length(suffix) == 1) time_stamp <- paste0(time_stamp, '_', suffix)

  return(time_stamp)
}
