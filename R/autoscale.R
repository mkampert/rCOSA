#' Mean-center and scale the column of a matrix
#'
#' @param x matrix
#' @param xmi matrix indicating which values are missing in \code{x}
#' @param xmo value that represents missings in the returned matrix
#' @param robust logical indicating whether to use \code{IQR(column_x) / 1.35} or just the standard deviation (\code{sd}).
#'
#' @return
#' a matrix in which the columns are mean-centered and scaled.#'
#'

autoscale <- function(x, xmi = NA, xmo = NA, robust = F) {
  y <- matrix(nrow = nrow(x), ncol = ncol(x))
  if (!is.na(xmi)) x[x == xmi] <- NA
  for (i in 1:ncol(x)) {
    u <- x[!is.na(x[, i]), i]
    if (robust) {
      scl <- IQR(u) / 1.349
      if (scl > 0) u <- (u - median(u)) / scl
      if (scl == 0) u <- (u - mean(u)) / sqrt(var(u))
    }
    else {
      scl <- var(u)
      if (scl > 0) {
        scl <- sqrt(scl)
        u <- (u - mean(u)) / scl
      }
    }
    y[!is.na(x[, i]), i] <- u
    y[is.na(x[, i]), i] <- xmo
  }
  invisible(y)
}
