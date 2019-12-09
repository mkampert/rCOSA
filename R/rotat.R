#' Rotate distance matrix (internal function rCOSA)
#'
#' @param xx distance matrix
#'
#' @return
#' Rotated matrix in which the first column has the largest eigen value etc.
#'

rotat <- function(xx) {
  rot <- eigen((t(xx)) %*% xx)
  rot <- rot[[2]]
  xx <- xx %*% rot
  return(xx)
}
