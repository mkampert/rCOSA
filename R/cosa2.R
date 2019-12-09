#' COSA 2 Dissimilarities
#'
#'  This function outputs a dissimilarity matrix of dissimilarities between the rows a data matrix computed by the COSA 2 algorithm. It is assumed that users are familiar with the COSA paper(s) or the vignette that comes with the rCOSA package, see references below.
#'
#'
#' @param X input data.frame, or matrix object in numeric mode. COSA calculates the dissimilarities for the rows in X.
#'
#' @param lX either an integer or a vector with as much elements as columns in \code{X}. This argument is ignored if \code{X} is a data.frame \cr
#' \cr.
#' Each element \code{lX[k]} represents a label for \code{X[,k]} as indicated below: \cr
#' 0 => ignore \code{X[, k]}, \cr
#' \cr
#' \code{1} => \code{X[, k]} is a numeric attribute, no target value, \cr
#' \code{2} => \code{X[, k]} is a numeric attribute, single target value, \cr
#' \code{3} => \code{X[, k]} is a numeric attribute, dual target values, \cr
#' \code{4} => \code{X[, k]} is a categorical attribute, no target value, \cr
#' \code{5} => \code{X[, k]} is a categorical attribute, single target value, \cr
#' \code{6} => \code{X[, k]} is a categorical attribute, dual target values. \cr
#' \cr
#' If only an integer scalar is give for \code{lX} than it is assumed that the integer value is the label for all columns in \code{X}. \cr Target values are specified in the targ, targ2 or qntls arguments (see below)\cr
#'
#' @param targ target values for computing targeted dissimilarities. The \code{targ[k]} is the target value for \code{X[,k]}. The value is ignored if \code{lX[k] = 0}, \code{1}, or \code{4}. \cr
#'         If \code{lX[k]=2} or \code{5}, then \code{targ[k]} contains the single target value. \cr
#'         If \code{lX[k]=3} or \code{6}, then \code{targ[k]} contains one of the two target values. \cr
#' \cr
#'         Special values: \cr
#'
#'            \code{targ = "low"} => use the \code{qntls[1]} quantile of \code{X[,k]} as a single target value. \cr
#'            \code{targ = "high"} => use the \code{qntls[2]} quantile of \code{X[,k]} as a single target value.    \cr
#'            \code{targ = "high/low"} => use the \code{qntls[1]} and \code{qntls[2]} quantiles of \code{X[,k]} as dual target values (see \code{qntls} argument described below). \cr
#'
#' @param targ2 the second target value when computing dual targeted dissimilarities. The \code{targ2[j]}  = second target value for \code{X[,k]}. The value is ignored if \code{lX[k] = 0},\code{1},\code{2},\code{4},\code{5}, or when \code{targ = "low"},\code{"high"}, or \code{"high/low"}.
#' \cr
#' @param knear size of number of objects in the near-neighborhoods which is used to calculate attribute weights for each object. By default \code{knear = sqrt(nrow(X))}, which, inside the function, is truncated to an  integer.
#' \cr
#' @param xmiss numeric value for missings in the data, by default it is set to \code{NULL}. In case you have coded missings in the data differently by using a number then indicate that number in this argument. The coded value for missings must be larger than any data value on any input variable.
#' \cr
#' @param lambda multiple attribute clustering incentive parameter. By default, the regularization parameter lambda is set to equal \code{0.2}.
#' \cr
#' @param qntls quantiles used for calculating high and/or low targets (ignored if \code{lX[k] = 0},\code{1},\code{4} or \code{targ} is NOT set to \code{"low"}, \code{"high"}, or \code{"high/low"}). \cr
#' \cr
#' \code{qntls[1]} = data quantile on each attribute used for low target \cr
#' \code{qntls[2]} = data quantile on each attribute used for high target \cr
#' \cr
#' NOTE: The \code{lX} vector and/or the value of \code{xmiss} can be specified as attributes to the input data matrix before invoking cosa \cr
#' \code{attr(X,"lX")<-lX} \cr
#' \code{attr(X,"xmiss")<-xmiss} \cr
#' \cr
#' When present the attribute values will be used whenever the corresponding arguments are missing. Specifying an argument value overrides the corresponding attribute values if present. If the attribute and its corresponding argument are both missing, the default values above are used.
#' \cr
#'
#' @param wtcomb by default is set to 'each', meaning that the maximum of the weights of object \code{X[i,]} and object \code{X[j,]} is choosen for each attribute \code{X[,k]}. One can also choose for the option 'whole', meaning that a whole weight vector is choosen of either object \code{i} or \code{j}, depending on which of the two vectors give the maximum weighted dissimilarity between object \code{i} and \code{j}.
#' \cr
#' @param relax the number with which the homotopy parameter eta should be incremented at each outer iteration (for more info see noit)
#' \cr
#' @param conv the convergence treshold that can reduce the maximum number of inner iterations.
#' \cr
#' @param niter the maximum number of inner iterations to stabilize the weights and dissimilarties given the homotopy parameter
#' \cr
#' @param noit the number of outer iterations ( make sure relax > 0 ) to transfer from the inverse exponential distance more closely to the sum of the weighted dissimilarities, obtained when a large enough number of outer iterations is chosen. Starting with the initial value of the homotopy parameter (equal to lambda) and using increments determined by relax one can calculate at what value the homotopy parameter will end.
#' \cr
#' @param stand equals \code{0} for no standardisation, \code{1} for robust standard scaling of the data, and \code{2} for standard scaling of the data. The defeault equals \code{1}.
#' \cr
#' @param pwr \code{1} for L_1 norm attribute distances, and \code{2} for L_2 norm attribute distances. The default equals \code{1}. This argument is ignored for categorical/ordinal attributes.
#'
#' @return
#' This function outputs a list that has as the first element the call, as the second element the dissimilarity matrix \code{out$D} of class dist, and by default also the weights of class matrix, and, if \code{crit} is set to \code{TRUE}, the values of the tuning parameters are also given in the output.
#'
#'
#' @references
#' Friedman, J. H. and Meulman, J. J. (2004). Clustering objects on subsets of attributes. \cr URL: \url{http://www-stat.stanford.edu/~jhf/ftp/cosa.pdf}
#' Kampert, M.M., Meulman J.J., Friedman J.H. (2017). rCOSA: A Software Package for Clustering Objects on Subsets of Attributes \cr URL: \url{https://link.springer.com/article/10.1007/s00357-017-9240-z}
#'
#' @author Maarten M.D. Kampert, Jacqueline J. Meulman, and Jerome H. Friedman. \cr Correspondence: \email{mkampert@@math.leidenuniv.nl}
#'
#'
#' @examples
#' data(ApoE3)
#' rslt_dflt_cosa2 <- cosa2(X = ApoE3)
#' # The weight of object 1 on attribute 1 in the NxP weight matrix W
#' rslt_dflt_cosa2$W[1, 1]
#' # COSA procedure for dual targeted dissimilarities:
#' rslt_duotrg_cosa <- cosa2(X = iris[, 1:4], lX = rep(3, ncol(iris[, 1:4])), targ = 'high/low')
#'
#'
#' @note
#' The output dissimilarity matrix can be used as input to most dissimilarity based clustering procedures in R in the same manner as the output of the procedure \code{\link[stats]{dist}}.
#'
#'
#' @seealso \code{\link[rCOSA]{hierclust}}, \code{\link[rCOSA]{getclust}}, and \code{\link[rCOSA]{smacof}}.
#'
#'
#' @export

cosa2 <- function(X, lX = NULL, targ = NULL, targ2 = NULL, knear = sqrt(nrow(X)),
                  xmiss = NULL, lambda = 0.2, qntls = c(0.05, 0.95), wtcomb = "each",
                  relax = 0.1, conv = 1e-05, niter = 1, noit = 100, stand = 1, pwr = 1) {

  OUT <- vector(length = 4, mode = "list")
  names(OUT) <- c("call", "D", "W", "tunpar")
  OUT$CALL <- sys.call()

  # Dimensions X:
  nrx <- as.integer(nrow(X))
  ncx <- as.integer(ncol(X))

  # Arrange the Targeting:
  if (is.null(targ)) {
    ktarg <- 1L
  } else if (is.numeric(targ)) {
    ktarg <- 0L
  } else if (targ == "high") {
    ktarg <- 1L
  } else if (targ == "high/low") {
    ktarg <- 1L
  } else if (targ == "low") {
    ktarg <- 2L
  } else {
    warning(paste(as.character(targ), "   invalid value for targ."))
  }

  if (is.null(targ2)) {
    ltarg <- 0L
  } else {
    ltarg <- 1L
  }

  tg <- cbind(
    if (ktarg == 0) {
      targ
    } else {
      rep(qntls[2], ncx)
    },
    if (ltarg) {
      targ2
    } else {
      rep(qntls[1], ncx)
    }
  )
  storage.mode(tg) <- "double"


  # Setting the Attribute Flags (1 - 6):

  seltarg <- if (!is.null(targ)) {
    2 * targ %in% c("low", "high") + 3 * (targ == "high/low")
  } else {
    1
  }

  if (length(lX) == ncx) {
    lXt <- lX
  } else if (length(lX) == 1) {
    lXt <- rep(lX, ncx)
    lX <- lXt
  } else if (is.null(lX) && (class(X) == "data.frame")) {
    lXt <- sapply(X, function(x) {
      if (class(x) %in% c("numeric", "double", "integer")) {
        out <- (1:3)[seltarg]
      } else if (class(x) %in% c(
        "factor", "character",
        "logical"
      )) {
        out <- (4:6)[seltarg]
      } else {
        stop("columns in X can only be of class factor or numeric")
      }
    })
    X <- sapply(X, function(x) as.numeric(x))
    lX <- lXt
  } else {
    stop("argument lX cannot be NULL when X is not a data.frame")
  }
  lX <- as.integer(lX)


  # Preparing the Data for Fortran (including missings):

  if (class(X) == "data.frame") X <- sapply(X, function(x) as.numeric(x))

  if (missing(xmiss) && !is.null(attr(X, "xmiss"))) {
    xmisst <- attr(X, "xmiss")
    X[is.na(X)] <- xmisst
  } else if (is.null(xmiss)) {
    xmisst <- max(9.9e+30, max(X, na.rm = TRUE) + 1)
    xmisst <- as.double(xmisst)
    X[is.na(X)] <- xmisst
  } else {
    xmisst <- as.double(xmiss)
  }
  storage.mode(X) <- "double"

  # Weighting strategy for dissimis...
  if (wtcomb == "each") {
    wtcom <- 1L
  } else if (wtcomb == "whole") {
    wtcom <- 2L
  } else {
    wtcom <- 1L
    warning(paste(as.character(wtcomb), " - invalid value for wtcomb: wtcomb='each' used."))
  }

  # All remaining integer parameters
  ipams <- c(niter, noit, stand, wtcom, ktarg, ltarg, itot = 0)
  ipams <- as.integer(ipams)
  kn <- as.integer(trunc(knear))

  # All remaining double parameters:
  dpams <- c(xmisst, lambda, qntls, relax, conv, pwr,
             crit1 = 0.0, crit2 = 0.0, eps = 0.0
  )
  dpams <- as.double(dpams)
  wi <- matrix(0, nrow = nrx, ncol = ncx)
  storage.mode(wi) <- "double"

  double.n2 <- double(nrx * (nrx - 1) / 2)

  Z <- .Fortran("cosa",
                no = nrx, ni = ncx, kn = kn, x = X,
                ip = ipams, dp = dpams, lx = lX, tg = tg,
                d = double.n2, wi = wi, PACKAGE = "rCOSA"
  )

  dis <- Z$d
  attr(dis, "Size") <- nrx
  class(dis) <- "dist"
  OUT$D <- as.dist(dis)
  # OUT$D = rCOSA:::normdist(OUT$D)

  OUT$W <- Z$wi
  # OUT$W = (ncx) * OUT$W/rowSums(OUT$W)

  OUT$tunpar <- c(
    Z$dp[8], Z$dp[2], Z$dp[10], Z$dp[9], kn,
    ipams[2], ipams[7], ipams[1], Z$dp[7], Z$dp[6], Z$dp[5]
  )
  names(OUT$tunpar) <- c(
    "crit", "lambda", "homotopy", "critrb",
    "Knn", "noit", "totit", "maxinit", "pwr", "conv", "relax"
  )

  # OUT$targets <- Z$tg

  return(OUT)
}
