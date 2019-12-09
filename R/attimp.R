#' Compute and plot relevance of the attributes for selected clusters
#'
#'    This function computes and displays the revelance (importance) of each attribute to the clustering of a selected group of objects. The relevance values are displayed in order of descending value. Optionally, on the same plot, are displayed corresponding ordered importances of randomly selected groups of objects of the same size.
#'
#'
#'
#' @param X data.frame or an input data matrix object of class matrix and mode numeric. Rows are objects and columns are attributes (variables).
#'
#' @param group vector of identities of objects in the group
#'
#' @param range range of ordered importance values plotted (all are computed)
#'
#' @param times number of randomly sampled groups for which ordered importances are also plotted
#'
#' @param lX either an integer or a vector of attribute flags with as much elements as columns in \code{X} (see the help file of \code{\link[rCOSA]{cosa2}}). This argument is ignored if \code{X} is a data.frame \cr
#' `
#' @param donames \code{TRUE} / \code{FALSE} => do/don't show names of attributes on plot.
#'
#' @param fast \code{TRUE} / \code{FALSE} => do/don't use fast approximate calculation (less stable estimate). Should only be done with very large groups (when regular calculation takes way too long).
#'
#' @param trans power used to transform important scale.
#'
#' @param maximp maximum obtainable importance value
#'
#' @param xmiss missing value flag (see the help file of \code{\link[rCOSA]{cosa2}}).
#'
#' @param horiz \code{TRUE}, \code{FALSE} => do/don't plot importances horizontally (ignored if time \code{> 0}).
#'
#' @param main plot title
#'
#' @param att.names vector of names for each attribute (ignored if donames = \code{FALSE}).
#'
#' @param ord.names a logical value that indicates whether the attributes labels should be ordered from high to low attribute importance (ignored if donames = \code{FALSE}).
#'
#' @param names.size expansion factor for names.
#'
#' @param ylim vertical range of plot.
#'
#' @param ntick number of tick marks on importance axis.
#'
#' @param rlines a logical value that indicates whether the replicate lines should show in the plot as well.
#'
#' @param targ target values for computing targeted dissimilarities. The \code{targ[k]} is the target value for \code{X[,k]}. The value is ignored if \code{lX[k] = 0}, \code{1}, or \code{4}. \cr
#'         If \code{lX[k]=2} or \code{5}, then \code{targ[k]} contains the single target value. \cr
#'         If \code{lX[k]=3} or \code{6}, then \code{targ[k]} contains one of the two target values. \cr
#' \cr
#'         Special values: \cr
#'
#'            \code{targ = "low"} \cr
#'            \code{targ = "high"} \cr
#'            \code{targ = "high/low"} \cr
#'
#' @param lwd line width (see \link[graphics]{par}).
#'
#' @param cache a logical value indicating whether the data should show as a variable (object) in the output as well (NB. NOT IN USE YET).
#'
#' @note
#' the \code{lX} vector and/or the value of \code{xmiss} can be specified as attributes to the input data matrix before invoking attimp (see \code{\link[rCOSA]{cosa2}})
#'
#' @return
#'    \item{$imp}{vector of all attribute importances ordered in descending value.}
#'    \item{$att}{vector of attribute identities in descending order of importance.}
#'
#' @examples
#' # make sure X and group are defined... For that, see getclust.
#' # at <- attimp(X, group[[2]])
#' # at <- attimp(X, group[[1]], 1:100, trans = 0.5, att.names = names, horiz = TRUE, main= 'Group 1')
#' # at <- attimp(X, group[[3]], 1:ncol(X), 10, main =' Group 3')
#' @seealso \code{\link[stats]{hclust}}, \code{\link[rCOSA]{getclust}}, and \code{\link[rCOSA]{attvalues}}.
#' @export

attimp <- function(X, group, range = 1:min(20, ncol(X)),
                   times = 0, lX = if (class(X) == "data.frame") {
                     lX <- NULL
                   } else {
                     lX <- 1
                   },
                   donames = FALSE, fast = FALSE, trans = 1, maximp = 20, xmiss = NULL, horiz = FALSE,
                   main = "Cluster", att.names = 1:ncol(X), ord.names = FALSE, names.size = 1,
                   ylim = c(0, wtave[o[1]]), ntick = 20, rlines = TRUE, targ = NULL, lwd = 1, cache = FALSE) {
  nrx <- nrow(X)
  ncx <- ncol(X)

  seltarg <- if (!is.null(targ)) {
    2 * targ %in% c("low", "high") + 3 * (targ == "high/low")
  } else {
    1
  }

  if (length(lX) == ncx) {
    lXt <- lX # labels X temporarily
  } else if (length(lX) == 1) {
    lXt <- rep(lX, ncx)
    lX <- lXt
  } else if (is.null(lX) && (class(X) == "data.frame")) {
    lXt <- sapply(X, function(x) {
      if (class(x) == "numeric") {
        out <- (1:3)[seltarg]
      } else if (class(x) %in% c("factor", "character", "logical")) {
        out <- (4:6)[seltarg]
      } else {
        stop("columns in X can only be of class factor or numeric")
      }
    })
    X <- sapply(X, function(x) as.numeric(x))
    lX <- lXt
    #    	rewrite X to matrix here!
  } else {
    stop("argument lX cannot be NULL when X is not a data.frame")
  }

  if (class(X) == "data.frame") X <- sapply(X, function(x) as.numeric(x))


  if (maximp <= 0) stop(" maximp must have a positive value.")
  if (missing(lX) && !is.null(attr(X, "lX"))) {
    lXt <- attr(X, "lX")
  } else {
    lXt <- lX
  }

  if (missing(xmiss) && !is.null(attr(X, "xmiss"))) {
    xmisst <- attr(X, "xmiss")
    X[is.na(X)] <- xmisst
  } else if (is.null(xmiss)) {
    xmisst <- max(9.9e+30, max(X, na.rm = TRUE) + 1)
    X[is.na(X)] <- xmisst
  } else {
    xmisst <- xmiss
  }

  eps <- 1 / maximp
  disp <- rep(0, ncol(X))
  wtave <- rep(0, ncol(X))
  for (k in 1:ncol(X)) {
    if (lXt[k] <= 0) {
      wtave[k] <- 0
      next
    }
    z <- X[X[, k] < xmisst, k]
    grp <- group[X[group, k] < xmisst]
    if (length(grp) < 3) {
      wtave[k] <- 0
      next
    }
    zg <- X[grp, k]
    if (lXt[k] < 4) {
      disp[k] <- IQR(z)
      if (disp[k] <= 0) {
        wtave[k] <- 0
        next
      }
      if (fast) {
        wtave[k] <- 1 / (IQR(zg) / disp[k] + eps)
      }
      else {
        wti <- rep(0, length(grp))
        for (i in 1:length(grp)) {
          wti[i] <- median(abs(zg[i] - zg))
        }
        wtave[k] <- 1 / (1.35 * mean(wti) / disp[k] + eps)
      }
    }
    else {
      u <- unique(z)
      s <- 0
      sq <- 0
      for (v in u) {
        t <- length(z[z == v])
        s <- s + t
        sq <- sq + t**2
      }
      if (s <= 0) {
        disp[k] <- 0
        wtave[k] <- 0
        next
      }
      disp[k] <- 1 - sq / s**2
      u <- unique(zg)
      s <- 0
      sq <- 0
      for (v in u) {
        t <- length(zg[zg == v])
        s <- s + t
        sq <- sq + t**2
      }
      if (s <= 0) {
        wtave[k] <- 0
      }
      else {
        wtave[k] <- 1 / ((1 - sq / s**2) / disp[k] + eps)
      }
    }
  }
  o <- order(-wtave)
  out <- list(imp = wtave[o], att = o)
  if (times <= 0) {
    if (horiz) {
      Xl <- "Importance"
      yl <- ""
    }
    else {
      Xl <- ""
      yl <- "Importance"
    }
    if (donames) {
      bplot(wtave[o[range]],
        names = as.character(att.names[if (ord.names) {
          o[range]
        } else {
          range
        }]),
        horiz = horiz, xlab = Xl, ylab = yl, zlim = ylim, trans = trans, ntk = ntick,
        names.size = names.size
      )
      title(main)
    }
    else {
      bplot(wtave[o[range]],
        ylab = yl, xlab = Xl, zlim = ylim,
        horiz = horiz, trans = trans, ntk = ntick
      )
      title(main)
      tickvalues <- c(0, seq(0, sqrt(length(range)), length.out = 6)[-1])^2
      axis(1, at = 1.2 * tickvalues, labels = tickvalues)
    }
  }
  else {
    plot(range, wtave[o[range]],
      ylab = "Importance", xlab = "Ordered Attributes", ylim = ylim,
      type = "n"
    )
    title(main)
    attimpline <- wtave[o[range]]
    wttave <- rep(0, length(wtave))
    for (it in 1:times) {
      mr <- trunc(nrow(X) * runif(length(group)) + 1)
      for (k in 1:ncol(X)) {
        if (lXt[k] <= 0) {
          wtave[k] <- 0
          next
        }
        if (disp[k] <= 0) {
          wtave[k] <- 0
          next
        }
        grp <- mr[X[mr, k] < xmisst]
        if (length(grp) < 3) {
          wtave[k] <- 0
          next
        }
        zg <- X[grp, k]
        if (lXt[k] < 4) {
          if (fast) {
            wtave[k] <- 1 / (IQR(zg) / disp[k] + eps)
          }
          else {
            wti <- rep(0, length(grp))
            for (i in 1:length(grp)) {
              wti[i] <- median(abs(zg[i] - zg))
            }
            wtave[k] <- 1 / (1.35 * mean(wti) / disp[k] + eps)
          }
        }
        else {
          u <- unique(zg)
          s <- 0
          sq <- 0
          for (v in u) {
            t <- length(zg[zg == v])
            s <- s + t
            sq <- sq + t**2
          }
          if (s <= 0) {
            wtave[k] <- 0
          }
          else {
            wtave[k] <- 1 / ((1 - sq / s**2) / disp[k] + eps)
          }
        }
      }
      wtave <- wtave[order(-wtave)]
      if (rlines) {
        lines(range, wtave[range], col = 3, lwd = lwd * 0.8, lty = 4)
      }
      wttave <- wttave + wtave
    }
    lines(range, wttave[range] / times, col = 2, lwd = lwd * 1.25)
    lines(range, attimpline, lty = 1, lwd = lwd * 1.5)
  }
  invisible(out)
}
