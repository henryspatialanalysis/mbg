#' Behavior when attaching the mbg package
#'
#' @details Yields a message if the INLA package namespace is not available.
#'
#' @param libname (character(1)) A character string giving the library directory where the
#'    package defining the namespace was found.
#' @param pkgname (character(1)) A character string giving the name of the package.
#'
#' @return (invisible) A message may be printed to the console.
.onAttach <- function(libname, pkgname) {
  # Check for INLA namespace
  if (!requireNamespace("INLA", quietly = TRUE)){
    msg <- paste0(
      "You are running with mbg package without INLA installed. You will not be able ",
      "to run the SPDE spatial models that rely on INLA functions. For details about ",
      "how to install INLA, see: https://www.r-inla.org/download-install"
    )
    packageStartupMessage(msg)
  }
}
