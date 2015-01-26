
#' COSA 2 Dissimilarities
#'
#'      *** ALPHA VERSION *** This function outputs a dissimilarity representation of the rows of a data matrix and corresponding weights computed by the COSA algorithm. It is assumed that users are familiar with the COSA paper(s) or the vignette that comes with the COSA package, see references below.
#'
#'
#' @param X input data.frame, or matrix object in numeric mode. COSA calculates the dissimilarities for the rows in X.
#'
#' @param lX either an integer or a vector with as much elements as columns in \code{X}. This argument is ignored if \code{X} is a data.frame \cr
#' \cr.
#' Each element \code{lX[k]} represents a label for \code{X[,k]} as indicated below: \cr
#' 0 => ignore \code{X[,k]}, \cr
#' \cr                                                                             
#' \code{1} => \code{X[,k]} is a numeric attribute, no target value, \cr
#' \code{2} => \code{X[,k]} is a numeric attribute, single target value, \cr  
#' \code{3} => \code{X[,k]} is a numeric attribute, dual target values, \cr
#' \code{4} => \code{X[,k]} is a categorical attribute, no target value, \cr
#' \code{5} => \code{X[,k]} is a categorical attribute, single target value, \cr      
#' \code{6} => \code{X[,k]} is a categorical attribute, dual target values. \cr
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
#' @param \code{targ2} the second target value when computing dual targeted dissimilarities. The \code{targ2[j]}  = second target value for \code{X[,k]}. The value is ignored if \code{lX[j] = 0},\code{1},\code{2},\code{4},\code{5}, or when \code{targ = "low"},\code{"high"}, or \code{"high/low"}. 
#' \cr
#' @param knear size (number of objects) of near-neighborhoods used to calculate attribute weights for each object. By default \code{knear = sqrt(nrow(X))}, which, in the function, is being truncated to an  integer.
#' \cr
#' @param xmiss numeric value for missings in the data, by default it is set to \code{NULL}. In case you have coded missings in the data differently by using a number then indicate that number in this argument. The coded value for missings must be larger than any data value on any input variable.
#' \cr
#' @param lambda multiple attribute clustering incentive parameter. By default, the regularization parameter lambda is set to equal \code{0.2}.
#' \cr
#' @param qntls quantiles used for calculating high and/or low targets (ignored if \code{lX[j] = 0},\code{1},\code{4} or \code{targ} is NOT set to \code{"low"}, \code{"high"}, or \code{"high/low"}). \cr
#' \cr
#' \code{qntls[1]} = data quantile on each attribure used for low target \cr
#' \code{qntls[2]} = data quantile on each attribure used for high target \cr
#' \cr                                                                          
#' NOTE: The \code{lX} vector and/or the value of \code{xmiss} can be specified as attributes to the input data matrix before invoking cosa \cr
#' \code{attr(X,"lX")<-lX} \cr
#' \code{attr(X,"xmiss")<-xmiss} \cr
#' \cr                          
#' When present the attribute values will be used whenever the corresponding arguments are missing. Specifing an argument value overrides the corresponding attribute values if present. If the attribute and its corresponding argument are both missing, the default values above are used.    
#' \cr
# @param wtcomb by default is set to 'each', meaning that the maximum of the weights of object \code{X[i,]} and object \code{X[j,]} is choosen for each attribute \code{X[,k]}. One can also choose for the option 'whole', meaning that a whole weight vector is choosen of either object \code{i} or \code{j}, depending on which of the two vectors give the maximum weighted dissimilarity between object \code{i} and \code{j}. <-- MK: "Needs a better description?"
#' \cr
# @param relax the number with which the homotopy parameter eta should be incremented at each outer iteration (for more info see noit)
#' \cr
# @param conv convergence treshold that can reduce the maximum number of inner iterations. <-- In a newer version of this package it will be specified how the number which is compared with this threshold is calculated.
#' \cr
#' @param niter the maximum number of inner iterations to stabilize the weights and dissimilarties given the homotopy parameter
#' \cr
#' @param noit the number of outer iterations (make sure relax > 0) to transfer from the inverse exponential distance more closely to the sum of the weighted dissimilarities, obtained when a large enough number of outer iterations is chosen. Starting with the initial value of the homotopy parameter (equal to lambda) and using increments determined by relax one can calculate at what value the homotopy parameter will end. <-- MK: "In a newer version this argument will hopefully be improved too."
#' \cr
#' @param stand equals `?\code{0} for no standardisation?, \code{1} for robust scaling of the data, and \code{2} for standard scaling of the data.
#' \cr
#' @param cosadir directory as a character string at which the executable COSA program can be found. By default it is specified as \code{NULL}, meaning that the executables in the COSA package are being used (which is advised).
#' \cr
#' @param wghts logical by default set to be \code{TRUE} meaning that the weights are given in the output
#' \cr
#' @param crit logical by default set to be \code{TRUE} meaning that the criterion and tuning parameters are given in the output.
#'
#' 
#' @return
#' This function outputs a list that has as the first element the call, as the second element the dissimilarity matrix \code{out$D} of class dist, and by default also the weights of class matrix, and, if \code{crit} is set to \code{TRUE}, the values of the tuning parameters are also given in the output. <-- MK: "In a newer version perhaps some other nice tweaks?"
#'
#'
#' @references
#' Friedman, J. H. and Meulman, J. J. (2004). Clustering objects on subsets of attributes. \cr URL: \url{http://www-stat.stanford.edu/~jhf/ftp/cosa.pdf}
#'
#'
#' @author Maarten M.D. Kampert, Jacqueline J. Meulman, and Jerome H. Friedman. \cr Correspondence: \email{mkampert@@math.leidenuniv.nl}
#'
#'
#' @examples
#' # Using the iris data as example:
#' out <- cosa2(X = iris[,1:4])
#' # Make a plot hierarchical clustering plot of the dissimilarities
#' hierclust(out$D)
#' # The weight of object 1 on attribute 1 in the NxP weight matrix W
#' out$W[1,1]
#' # COSA procedure for dual targeted dissimilarities + transformed data is being asked using sdat
#' out <- cosa2(X = iris[,1:4],lX = rep(3,ncol(iris[,1:4])),targ='high/low')                       
#'
#'
#' @note
#' The output dissimilarity matrix can be used as input to most dissimilarity based clustering procedures in R in the same manner as the output of the procedure \code{\link[stats]{dist}}.
#'
#'
#' @seealso \code{\link[stats]{hclust}} and other functions that will be mentioned during the development of this package.
#'
#'
#' @export

cosa2 <- function (X, lX = if(class(X) == 'data.frame'){NULL} else {1}, targ = NULL, targ2 = NULL, knear = sqrt(nrow(X)), 
    xmiss = NULL, lambda = 0.2, qntls = c(0.05, 0.95), wtcomb = "each", 
    relax = 0.1, conv = 1e-05, niter = 5, noit = 100, stand = 1, 
    cosadir = NULL, wghts = TRUE, crit = TRUE, clean = TRUE, fnameprefix = NULL, pwr = 1, ...) 
{

#	OUT <- list(call = match.call())
    OUT <- list(call = sys.call())	


	# Find the executable and OS platform
    OSx = grep(c("x86_64"), R.version$platform) && grep(c("apple"), 
        R.version$platform)
    linux = grep(c("x86_64"), R.version$platform) && grep(c("linux"), 
        R.version$platform)
    win32 = grep(c("i386"), R.version$platform) && grep(c("mingw32"), 
        R.version$platform)
    win64 = grep(c("x86_64"), R.version$platform) && grep(c("mingw32"), 
        R.version$platform)
    platform = c("OSx", "linux", "win32", "win64")[!c(is.na(OSx), 
        is.na(linux), is.na(win32), is.na(win64))]
    if (is.null(cosadir)) {
        cosadir = paste(searchpaths()[search() == "package:rCOSA"], 
            "/execs/", platform, sep = "")
    }


 
    # Targeting section
    if (is.null(targ)){
        ktarg = 1
    } else if (is.numeric(targ)) {
        ktarg = 0
    } else if (targ == "high") {
        ktarg = 1
    } else if (targ == "high/low") {
        ktarg = 1
    } else if (targ == "low") {
        ktarg = 2
    } else {
        warning(paste(as.character(targ), "   invalid value for targ."))
    }
    
    if (is.null(targ2)) {
        ltarg = 0
    } else {
        ltarg = 1
    }


	# Variables labels and Data section

    nrx <- nrow(X)
    ncx <- ncol(X)

    seltarg <- if(!is.null(targ)){
        2*targ%in%c('low','high') + 3*(targ =='high/low')
    	} else {
    		1
    	}
    	
    if(length(lX) == ncx){
        lXt <- lX     # labels X temporarily    	
    } else if(length(lX) == 1) {
        lXt <- rep(lX, ncx)
        lX <- lXt
    } else if (is.null(lX) && (class(X) == 'data.frame')){
        lXt <- sapply(X, function(x){
        	if(class(x) %in% c('numeric','double','integer')){
        			out <- (1:3)[seltarg]
        		} else if (class(x) %in% c('factor', 'character', 'logical')){
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
    
    if(class(X)=='data.frame') X <- sapply(X, function(x) as.numeric(x))

    # Missing section
    
    if (missing(xmiss) && !is.null(attr(X, "xmiss"))) {
        xmisst = attr(X, "xmiss")
        X[is.na(X)] <- xmisst
    } else if (is.null(xmiss)){
    	xmisst <- max(9.9e+30, max(X, na.rm = TRUE) + 1)
        X[is.na(X)] <- xmisst
    } else {
    	xmisst = xmiss
    }
    
    
    # weight combination section    
    
    if (wtcomb == "each") {
        wtcom <- 1
    } else if (wtcomb == "whole") {
        wtcom <- 2
    } else {
        wtcom <- 1
        warning(paste(as.character(wtcomb), " - invalid value for wtcomb: wtcomb='each' used."))
    }
    
    
    ic <- c(nrx, ncx, trunc(knear), ktarg, niter, wtcom, ltarg, 
        noit, stand)


    	zz <- file(paste(cosadir, "/",fnameprefix,"iparms.cos", sep = ""), "wb")
	    writeBin(as.integer(ic), zz, size = 4)
	    close(zz)

    	zz <- file(paste(cosadir, "/",fnameprefix,"lx.cos", sep = ""), "wb")
	    writeBin(as.integer(lXt), zz, size = 4)
    	close(zz)

	    zz <- file(paste(cosadir, "/",fnameprefix,"rparms.cos", sep = ""), "wb")
    	writeBin(as.double(c(xmisst, lambda, qntls, relax, conv, pwr)), 
        zz, size = 4)
	    close(zz)

    	zz <- file(paste(cosadir, "/",fnameprefix,"data.cos", sep = ""), "wb")
	    writeBin(as.double(X), zz, size = 4)
    	close(zz)

	    if (ktarg == 0) {
    	    zz <- file(paste(cosadir, "/",fnameprefix,"targs.cos", sep = ""), "wb")
        	writeBin(as.double(targ), zz, size = 4)
	        close(zz)
    	}
	    if (ltarg == 1) {
    	    zz <- file(paste(cosadir, "/",fnameprefix,"dtargs.cos", sep = ""), "wb")
        	writeBin(as.double(targ2), zz, size = 4)
	        close(zz)
    	}
	    wd <- getwd()
    	
    	
    	
    	setwd(cosadir)

		  #CLI usage: ./COSA_2 iparms.cos rparms.cos lx.cos data.cos targs.cos dtargs.cos dist.cos weights.cos tunpar.cos

		inoutfiles <- c("iparms.cos", "rparms.cos", "lx.cos", "data.cos", "targs.cos", "dtargs.cos", "dist.cos", "weights.cos", "tunpar.cos") 
		extra_args <- paste(fnameprefix, inoutfiles, sep = "", collapse = " ") 

	    switch(platform, 
	    	OSx = system(paste("./COSA_2", extra_args, sep = " "), ...),
	    	linux = system(paste("./COSA_2", extra_args, sep = " "), ...), 
    		win32 = system(paste("cmd.exe /c COSA_2.exe", extra_args, sep = " "), minimized = F, invisible = FALSE, show.output.on.console = FALSE, ...), 
        win64 = system(paste("cmd.exe /c COSA_2.exe", extra_args, sep = " "), minimized = F, invisible = FALSE, show.output.on.console = FALSE, ...)
        
        )

	    setwd(wd)
	    	    
    	zz = file(paste(cosadir, "/",fnameprefix,"dist.cos", sep = ""), "rb")
	    dis = readBin(zz, numeric(), size = 4, n = nrx * (nrx - 1)/2)
    	close(zz)
	    attr(dis, "Size") = nrx
    	class(dis) <- "dist"
	    OUT$D = as.dist(dis)
    	if (wghts) {
        	zz = file(paste(cosadir, "/",fnameprefix,"weights.cos", sep = ""), "rb")
	        weights = readBin(zz, numeric(), size = 4, n = nrx * ncx)
	        close(zz)
    	    weights = matrix(weights, ncol = ncx, nrow = nrx)
        	weights = weights * (ncx)
	        OUT$W = weights
    	}
	    #if (sdat) {
    	#    zz = file(paste(cosadir, "/",fnameprefix,"sdata.cos", sep = ""), "rb")
        #	sdata = readBin(zz, numeric(), size = 4, n = nrx * ncx)
	    #    close(zz)
    	#    DAT = matrix(sdata, ncol = ncx, nrow = nrx)
        # 	DAT[DAT == xmiss] = NA
	    #     OUT$DAT = DAT
   	    #}
    
	    if (crit) {
    	    zz = file(paste(cosadir, "/",fnameprefix,"tunpar.cos", sep = ""), "rb")
        	tunpar = readBin(zz, numeric(), size = 4, n = 7)
	        close(zz)
    	    #tunpar[5:7] = as.integer(tunpar[5:7]) # has no use.... all elements should have the same class.
	        names(tunpar) = c("crit", "lambda", 'homotopy', 'RMSD', 'Knn', 'noit', 'totit')
    	    #print(tunpar)
	        OUT$tunpar = as.list(tunpar)
    	}
    	
 	if (clean) {
	    system(paste("rm ", cosadir, "/",fnameprefix,"*.cos", sep = ""))
	    #system(paste("rm ", cosadir, "/",fnameprefix,"*.666", sep = ""))
    }
    return(OUT)
}



