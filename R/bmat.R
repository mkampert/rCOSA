#' Create the so-called B matrix
#'
#' An internal function int the \link[rCOSA]{smacof} algorithm.
#'
#' @param odist matrix of observed distances
#' @param fdist matrix of fitted distances
#'
#' @return
#' the B matrix
#'
bmat <- function(odist, fdist) {
  n <- nrow(dist)
  nnn <- matrix(0, nrow = n, ncol = n)
  bb <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    nnn[i, i] <- 1
  }
  for (i in 1:n)
  {
    sb <- 0
    for (j in 1:n)
    {
      if (fdist[i, j] != 0) {
        bb[i, j] <- odist[i, j] / fdist[i, j]
        sb <- sb - bb[i, j]
      }
    }
    bb[i, i] <- bb[i, i] + sb
  }
  bmat <- -bb
  return(bmat)
}
