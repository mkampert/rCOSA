#' SMACOF
#'
#' A metric multidimensional scaling algorithm that minimizes a least squares loss function based on the dissimilarities Heiser \& De Leeuw, 80 and 77; Gutman 68). The lettters in SMACOF stand for Scaling by MAjorizing a COmplicated Function.
#' 
#' @param D dissimilarities/distances of class dist.
#' @param interc logical for including an additive constant in the smacof algorithm.
#' @param inicon intial configuration. By default this is specified as NULL, here the SMACOF algorithm specifies...
#' @param niter max number of iterations to use, default is set to 20: 'to iterate is heaven, to converge is divine'. Convergence does not necessarily give you the best non-linear mapping of the dissimilarity structure.
#' @param groupnr groupnr when known, each object can be given a number to which group it belongs
#' @param colv vector containing color names for each group number
#' @param main title of the plot, by default there is no title. 
#' @param pchset plotting 'character', i.e. symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols. The full set of S symbols is available with \code{pch = 0:18}. For more information see the help file of the function \code{\link[graphics]{points}}.
#' @param ... arguments which can be parsed to either the \code{\link[mva]{hclust}} or \code{\link[stats]{cmdscale}} procedures.
#' @param PLOT whether to plot or not to plot, by default TRUE.
#' @param VERBOSE whether to give a stdout on the criterion for each iteration, by default set to equal TRUE.
#' @return 
#' Similar output as \code{\link[stats]{cmdscale}}
#'
#' @export                                                                               
                                        
smacof <- function(D,niter = 100, interc = 1, inicon = NULL, groupnr = NULL, colv = palette()[c(8,2,4,3,5,6,7,1)], main = "Multidimensional Scaling", k = 2, pch = 16, PLOT = TRUE, VERBOSE = TRUE, ...)
{
        dhat<-ddelta<-D<-normdist(D)
        n<-dimdata(D)
        delta<-lowmat(D)
        ycol<-matrix(c(0,0),nrow=1,byrow=T)

	    if (is.null(inicon)) {
        	inicon <- cmdscale(D, k)
    	}
	    yy <- xx <- inicon

        dd<-lowmat(dist(xx), ...)
        str0 <- str <-sum((delta-dd)**2)
        for (i in 1:niter)
        {
        if (niter !=0)
           {
            xx<-(1/n)*(bmat(delta,dd))%*%xx
            dd<-dist(xx)
            if (interc==1)       
               {
               ddd<-dd
               dhat<-lm(dd~ddelta)$fitted.values;delta<-lowmat(normdist(dhat))
               if ((a<-min(dhat))<=0)
                  {
                 ddd<-.01+dd+abs(a);dhat<-lm(ddd~ddelta)$fitted.values;delta<-lowmat(normdist(dhat))
                  }
               }
        dd<-lowmat(dist(xx))
            #str=(sum((delta-dd)**2))**.5
            str <- snorm <-fit(delta,dd)
	        if(VERBOSE) cat("it = ", i, ",  crit = ", str, " \n")
	        if( abs(str0 - str) < 1e-6 ) break
	        str0 <- str
            }
    }
        xx <- rotat(xx)
        
	    if(!is.null(groupnr)){
	    		cols <- colv[groupnr + 1]
	    		} else {
	        cols <- 'black'	    			
	    		}

        if (PLOT){
            lims <- c( 
            ll = floor(100 * min( c(xx[, 1], xx[, 2]) ) )/100,        	
            ul = ceiling(100 * max(c(xx[, 1], xx[, 2]) ) )/100
          )  
          
	        	 # plot(xx,asp=1,type='n',xlab='Dimension 1 ',ylab='Dimension 2')
    	     # points(xx,asp=1,col="gray",pch=21,bg="gray")
    	     plot(xx, xlab = "Dimension 1", ylab = "Dimension 2",
    	    	main = main, xlim = lims, ylim = lims,
    	    	xaxs = "i", yaxs = "i", pch = pch, col = cols, ...)
        }

    	delta<- normdist(dhat)
    	return(list(D=delta, X=xx, interc = interc))	
    	
}




# de Leeuw J, Heiser WJ (1977). ?Convergence of Correction Matrix Algorithms for Multidi- mensional Scaling.? In J Lingoes (ed.), Geometric Representations of Relational Data, pp. 735?752. Mathesis Press, Ann Arbor.

# de Leeuw J, Heiser WJ (1980). ?Multidimensional Scaling with Restrictions on the Configura- tion.? In P Krishnaiah (ed.), Multivariate Analysis, Volume V, pp. 501?522. North Holland Publishing Company, Amsterdam.

# Guttman L (1968). ?A General Nonmetric Technique for Fitting the Smallest Coordinate Space for a Configuration of Points.? Psychometrika, 33, 469?506.

# literature on additive constant?
