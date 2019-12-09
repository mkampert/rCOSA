#' Return matrix with only the lower diagonal
#'
#' @param udata a distance object, or vector with distances
#'
#' @return
#' matrix with only the lower diagonal
#'

lowmat <- function(udata) {
  n <- dimdata(udata)
  lowmat <- matrix(0, nrow = n, ncol = n)
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in i:(n - 1)) {
      k <- k + 1
      lowmat[i, j + 1] <- udata[k]
      lowmat[j + 1, i] <- udata[k]
    }
  }
  return(lowmat)
}
