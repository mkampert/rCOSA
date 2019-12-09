#' Compute fit between two distance matrices
#'
#' @param delta observed distance matrix
#' @param dd euclidean distance matrix of a configuration
#'
#' @return
#' the value for the fit (STRESS)
#'
#'
fit <- function(delta, dd) {
  stot <- sum(delta**2)
  craw <- sum(delta * dd)
  traw <- sum(dd**2)
  sraw <- stot + traw - 2 * craw
  snor <- (1 - ((craw**2) / (traw * stot)))**.5
  #    print(c(stot,craw,traw,sraw,snor))
  return(snor)
}
