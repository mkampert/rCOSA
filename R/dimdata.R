#' Retrieve the number of observations from a distance object
#'
#' @param data an object of class \code{dist}, or a vector with distances
#'
#' @return
#' the number of observation in the data set

dimdata <- function(data) {
  dimdata <- 1 + sqrt(2 * length(data) + .25) - .5
  return(dimdata)
}
