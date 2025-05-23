% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logging.R
\name{logging_start_timer}
\alias{logging_start_timer}
\title{Start logging timer}
\usage{
logging_start_timer(msg, echo = TRUE, indentation_text = "  ")
}
\arguments{
\item{msg}{(\code{character(1)}) Logging message}

\item{echo}{(\code{logical(1)}, default TRUE) Should the message be written to screen?}

\item{indentation_text}{(\code{character(1)}, default "  ") Text that will be repeated at
the beginning of the message for each layer of indentation}
}
\description{
Start a nested timer with an optional message
}
\examples{
mbg::logging_start_timer(msg = 'Test logging')
Sys.sleep(0.1)
mbg::logging_stop_timer()
log_results <- mbg::logging_get_timer_log()
print(log_results)

}
