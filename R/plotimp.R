#' Replot attribute relevances with different plot parameters (after calling "attimp")
#'
#' This function allows replotting of the attribute revelances (importances) with different display parameters after calling "attimp", without having to recompute the relevance values. It also allows flexibility in how the results are displayed.
#'
#' @param attimp output from "attimp"
#'
#'
#' @param atts vector of attribute identities. This vector can be numeric, referencing the respective attribute numbers, or character, referencing the respective attribute names. If the latter, then the parameter "att.names" (see below) must be supplied.
#'
#' @param imporder \code{TRUE} / \code{FALSE}  => do/don't plot relevance values in order of descending value.
#' @param donames \code{TRUE} / \code{FALSE} => do/don't show names of attributes on plot.
#' @param trans power used to transform importance scale.
#' @param horiz \code{TRUE} / \code{FALSE} => do/don't plot importances horizontally
#' @param main plot title
#' @param att.names vector of names for each attribute (ignored if donames= \code{FALSE}).
#' @param names.size expansion factor for names.
#' @param plot \code{TRUE} / \code{FALSE}  => do/don't display barplot.
#' @param ylim vertical range of plot.
#' @param ntick number of tick marks on importance axis.
#'
#'
#'
#' @return
#' values of attribute importances in requested order.
#'
#' @examples
#' # make sure an at object is available...
#' # plotimp(at)
#' # plotimp(at,1:100,trans=0.5,horiz=TRUE)
#' # plotimp(at,c(27,43,10,16))
#' # plotimp(at,c('age','income','occupation'),att.names=demographics)
#' @seealso
#' \code{\link[rCOSA]{attimp}}
#'
#'
#' @export
plotimp <-
  function(attimp, atts = attimp$att[1:min(40, length(attimp$att))],
           imporder = F, donames = T, trans = 1, horiz = F, main = "Cluster",
           att.names = 1:length(attimp$att), names.size = 1, plot = T,
           ylim = c(0, max(imps)), ntick = 20) {
    if (is.character(atts)) {
      if (missing(att.names)) {
        stop(' addressing attributes by name requires input "att.names"')
      }
      attnums <- match(atts, att.names)
    }
    else {
      attnums <- atts
    }
    imps <- attimp$imp[match(attnums, attimp$att)]
    o <- 1:length(atts)
    if (imporder) o <- order(-imps)
    if (plot) {
      if (horiz) {
        xl <- "Importance"
        yl <- ""
      }
      else {
        xl <- ""
        yl <- "Importance"
      }
      if (donames) {
        bplot(imps[o],
          names = as.character(att.names[attnums[o]]),
          horiz = horiz, xlab = xl, ylab = yl, zlim = ylim, trans = trans, ntk = ntick,
          names.size = names.size
        )
        title(main)
      } else {
        bplot(imps[o],
          ylab = yl, xlab = xl, zlim = ylim,
          horiz = horiz, trans = trans, ntk = ntick
        )
        title(main)
      }
    }
    if (plot) {
      invisible(imps)
    } else {
      imps
    }
  }
