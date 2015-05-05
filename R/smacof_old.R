#' SMACOF
#'
#' A metric multidimensional scaling algorithm that minimizes a least squares loss function based on the dissimilarities Heiser \& De Leeuw, 80 and 77; Gutman 68). The lettters in SMACOF stand for Scaling by MAjorizing a COmplicated Function.
#' 

#' @param D dissimilarities/distances of class dist.
#' @param norm whether to normalize the dissimilarties. The default is set equal to \code{TRUE}
#' @param inicon intial configuration. By default this is specified as NULL, here the SMACOF algorithm specifies...
#' @param niter max number of iterations to use, default is set to 20: 'to iterate is heaven, to converge is divine'. Convergence does not necessarily give you the best non-linear mapping of the dissimilarity structure.
#' @param groupnr groupnr when known, each object can be given a number to which group it belongs
#' @param colv vector containing color names for each group number
#' @param main title of the plot, by default there is no title. 
#' @param pchset plotting 'character', i.e. symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols. The full set of S symbols is available with \code{pch = 0:18}. For more information see the help file of the function \code{\link[graphics]{points}}.
#' @param ... arguments which can be parsed to either the \code{\link[mva]{hclust}} and \code{\link[stats]{cmdscale}} procedures.
#' @param PLOT whether to plot or not to plot, by default TRUE.
#' @param VERBOSE whether to give a messages on the criterion for each iteration, by default set to equal TRUE.
#' @return 
#' Similar output as \code{\link[stats]{cmdscale}}
#'
#' @export                                                                               
                                        
smacof_old <- function (D, normalize = TRUE, inicon = NULL, niter = 20, interc = TRUE,
    groupnr = NULL, colv = palette(), main = NULL, pchset = 20, 
    ..., PLOT = TRUE, VERBOSE = TRUE){
    n <- (1 + sqrt(2 * length(D) + 0.25) - 0.5)
    if (normalize) {
        fac <- sum(D^2)
        D <- D * ((n/fac)^0.5)
    }
    delta <- as.matrix(D)
    ycol <- matrix(c(0, 0), nrow = 1, byrow = TRUE)
    if (is.null(inicon)) {
        inicon <- cmdscale(D)
    }
    xx = inicon
    dd = as.matrix(dist(xx))
    crit <- sum((delta - dd)^2)
    for (i in 1:niter) {
        if (niter != 0) {
            bmat <- -as.matrix(delta/dd)
            diag(bmat) <- -rowSums(bmat, na.rm = TRUE)
            xx <- (1/n) * (bmat) %*% xx
            if (interc==1){
		        dhat<-lm(dd~delta)$fitted.values
		        if (normalize) {
     		 	   fac <- sum(dhat^2)
			        delta <-  dhat * ((n/fac)^0.5)
			        }
		    }
            dd <- as.matrix(dist(xx))
            crit <- (sum((delta - dd)^2))^0.5
        }
        xx <- xx %*% eigen((t(xx) %*% xx))[[2]]
        if(VERBOSE) cat("it = ", i, ",  crit = ", crit, " \n")
    }
        cols <- 'black'
	    cols <- if(!is.null(groupnr)) colv[groupnr + 1] 
    	lims <- c( 
    		ll = floor(100 * min( c(xx[, 1], xx[, 2]) ) )/100,        	
    	    ul = ceiling(100 * max(c(xx[, 1], xx[, 2]) ) )/100
    	    )  
        	
	if (PLOT){
	    plot(x = xx[, 1], y = xx[, 2], xlab = "Dimension 1", 
	    	ylab = "Dimension 2", main = main, col = cols, 
	    	xlim = lims, ylim = lims, xaxs = "i", 
	    	yaxs = "i", pch = pchset, ...)
	    }
	    invisible(list(xx, lims))

}

# de Leeuw J, Heiser WJ (1977). ?Convergence of Correction Matrix Algorithms for Multidi- mensional Scaling.? In J Lingoes (ed.), Geometric Representations of Relational Data, pp. 735?752. Mathesis Press, Ann Arbor.

# de Leeuw J, Heiser WJ (1980). ?Multidimensional Scaling with Restrictions on the Configura- tion.? In P Krishnaiah (ed.), Multivariate Analysis, Volume V, pp. 501?522. North Holland Publishing Company, Amsterdam.

# Guttman L (1968). ?A General Nonmetric Technique for Fitting the Smallest Coordinate Space for a Configuration of Points.? Psychometrika, 33, 469?506.
