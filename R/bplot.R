#' Barplot for the rCOSA results
#'
#'    Create a Barplot for rCOSA results
#'
#' @param x a vector or matrix of values
#' @param xlab label on the horizontal axis
#' @param ylab label on the vertical axis
#' @param zlim limits for the axis representing the height of the bars
#' @param names label on the vertical axis
#' @param horiz logical value that indicates whether the bars are shown horizontal or not (default is \code{FALSE}).
#' @param trans ...
#' @param ntk number of equally spaced thicks on the barplot axis
#' @param names.size  sizes of the names (expressed in values for \code{cex.names})

bplot <- function(x, xlab = "", ylab = "", zlim = c(0, 1), names = NULL, horiz = F, trans = 1, ntk = 20, names.size = 1) {
  if (trans == 1) {
    if (horiz) {
      barplot(x,
        names = names, horiz = T, xlab = xlab, ylab = ylab, xlim = zlim, las = 1,
        cex.names = names.size
      )
    }
    else {
      barplot(x,
        names = names, horiz = F, xlab = xlab, ylab = ylab, ylim = zlim, las = 2,
        cex.names = names.size
      )
    }
  }
  else if (horiz) {
    barplot(x**trans,
      names = names, horiz = T, xlab = xlab, ylab = ylab, las = 1,
      xlim = zlim**trans, axes = F, cex.names = names.size
    )
    lab <- pretty(x, n = ntk)
    axis(side = 1, at = lab**trans, labels = as.character(lab))
  }
  else {
    barplot(x**trans,
      names = names, xlab = xlab, ylab = ylab, ylim = zlim**trans, axes = F,
      las = 2, cex.names = names.size
    )
    lab <- pretty(x, n = ntk)
    axis(side = 2, at = lab**trans, labels = as.character(lab))
  }
}
