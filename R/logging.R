## #######################################################################################
##
## MBG package logging utilities
##
## #######################################################################################

#' Add global indent for logging
#' 
#' @description Change a global option that controls indentation for MBG package logging
#' 
#' @details Increases the `MbgLoggingIndent` option, with default of one if not set
#' 
logging_add_indent <- function(){
  options(MbgLoggingIndent = getOption('MbgLoggingIndent', default = 0) + 1)
  invisible()
}

#' Remove global indent for logging
#' 
#' @description Change a global option that controls indentation for MBG package logging
#' 
#' @details Decreases the `MbgLoggingIndent` option, with default of zero if not set
#' 
logging_drop_indent <- function(){
  options(MbgLoggingIndent = getOption('MbgLoggingIndent', default = 1) - 1)
  invisible()
}


#' Start logging timer
#' 
#' @description Start a nested timer with an optional message
#' 
#' @param msg (`character(1)`) Logging message
#' @param echo (`logical(1)`, default TRUE) Should the message be written to screen?
#' @param indentation_text (`character(1)`, default "  ") Text that will be repeated at
#'   the beginning of the message for each layer of indentation
#' 
#' @importFrom tictoc tic
#' @export
logging_start_timer <- function(msg, echo = TRUE, indentation_text = '  '){
  # Add indentation to the message
  indentation <- rep(indentation_text, times = getOption('MbgLoggingIndent', default = 0)) |>
    paste0(collapse = '')
  msg <- paste0(indentation, msg)
  tictoc::tic(msg = msg)
  if(echo) message(msg)
  # Add indentation for future nested timers
  logging_add_indent()
  invisible()
}


#' End logging timer
#' 
#' @description End a nested timer
#' 
#' @param echo (`logical(1)`, default = TRUE) Should the message be written to screen?
#' 
#' @importFrom tictoc toc
#' @export
logging_stop_timer <- function(echo = TRUE){
  tictoc::toc(log = TRUE, quiet = !echo)
  # Remove nested indentation from this timer
  logging_drop_indent()
  invisible()
}


#' Get timer log
#' 
#' @description Return a log of all timed events as a data.table
#' 
#' @param clear_log (`logical(1)`, default FALSE) Should the log be cleared afterwards?
#' @param deindent (`logical(1)`, default TRUE) Should leading whitespace be removed from
#'   timer messages?
#' 
#' @importFrom tictoc tic.log tic.clearlog
#' @import data.table
#' @export 
logging_get_timer_log <- function(clear_log = FALSE, deindent = TRUE){
  timer_log <- tictoc::tic.log(format = FALSE) |>
    lapply(data.table::as.data.table) |>
    data.table::rbindlist(fill = T, use.names = T)
  if('tic' %in% colnames(timer_log)) timer_log <- timer_log[order(tic)]
  if(all(c('tic', 'toc') %in% colnames(timer_log))) timer_log[, elapsed := toc - tic]
  if(deindent){
    if('msg' %in% colnames(timer_log)) timer_log[, msg := gsub('^ +', '', msg)]
    if('callback_msg' %in% colnames(timer_log)){
      timer_log[, callback_msg := gsub('^ +', '', callback_msg)]
    }
  }
  if(clear_log) tictoc::tic.clearlog()
  return(timer_log)
}
