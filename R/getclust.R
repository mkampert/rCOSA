



#' Identify Clusters in a Dendrogram
#'
#'    Reads the position of the graphics pointer when the left mouse button is pressed. It then cuts the tree at the vertical position of the pointer and draws a colored rectangle (rec.col) around the cluster containing the horizontal position of the pointer, and identifies the objects in the corresponding cluster. This can be repeated up to a maximum of N times. Each press of the left mouse button selects a (different) cluster. The new cluster is colored rec.col and the color of the previous one is changed to old.col. When finished, pressing the right mouse button and selecting 'stop' returns control to the command window.
#'
#'
#'
#'
#' @param hc an object of the type produced by \code{\link[rCOSA]{hierclust}} or \code{\link[stats]{hclust}}
#' @param ngr maximum number of clusters that can be selected
#' @param rec.col color of the most recently selected cluster
#' @param old.col color of the previously selected clusters.
#'
#'
#' @return
#' (Invisibly) returns a list where each element contains a vector of object identities contained in the respectively selected clusters.
#'
#' @examples
#' # need to have an object hc produced by hierclust() or hclust()
#' # groups <- getclust(hc, NrG= 3)
#' @note
#' Make sure you have opened a plot device of your dendrogram already! The vector of object identities for each selected cluster are contained in \code{group[[1]]}, \code{group[[2]]}, \code{...}, \code{group[[NrG]]} in order of selection, where \code{NrG} is the number of selected groups (clicks of the left mouse button). Press Esc or click the right mouse button to end the selection procedure.
#'
#'
#' @seealso \code{\link[rCOSA]{hierclust}}, \code{\link[stats]{hclust}}, and \code{\link[cluster]{pam}}
#'
#'
#' @export

getclust <- function(hc, ngr = 20, rec.col = "blue", old.col = "blue") {
  # plan: modify this function to make automatic clustering possible too. Using pam?
  retval <- list()
  oldk <- NULL
  oldx <- NULL
  cat("Showing dynamic visualisation. Press Escape/Ctrl + C to stop.")
  for (n in 1:ngr) {
    x <- locator(1)
    if (is.null(x)) {
      break
    }
    k <- min(which(rev(hc$height) < x$y))
    k <- max(k, 2)
    if (!is.null(oldx)) {
      rect.hclust(hc, k = oldk, x = oldx, border = old.col)
    }
    retval[[n]] <- unlist(rect.hclust(hc,
      k = k, x = x$x,
      border = rec.col
    ))
    oldx <- x$x
    oldk <- k
    ngr <- n
  }

  grps <- rep(0, length(hc$order))

  for (n in 1:ngr) {
    grps[retval[[n]]] <- n
  }
  names(retval) <- paste("grp", 1:n, sep = "")

  invisible(list(grps = grps, index = retval))
}
