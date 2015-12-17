#' Plot distributions of selected attributes within a specified cluster       
#'
#' This procedure plots a matrix of histograms. Each column pertains to one of the selected attributes. The first row shows the distribution for those objects within the specified group. The second row plots the corresponding distribution for those objects not in the specified group. Optionally (\code{all= TRUE}), a third row of histograms is produced showing the distribution of the corresponding attribute over all of the data. In all histograms, only nonmissing values of the respective attributes are plotted.
#'
#'
#' @param X input data. Either a data.frame or a matrix in mode numeric. Rows are objects and columns are attributes (variables).
#' @param group vector of identities of objects in the group.
#' @param atts vector of identities of selected attributes.
#' @param lX vector of attribute flags (see \code{\link[COSA]{cosa}}).
#' @param att.names vector of names for each attribute.
#' @param group.name name for selected group (cluster).
#' @param xmiss missing value flag (see \code{\link[COSA]{cosa}}).
#' @param nbins (approXimate) number of histogram bins.
#' @param all \code{TRUE}, \code{FALSE} => do/don't show histograms for all data (third row).
#' 
#' @note
#' The \code{lX} vector and/or the value of \code{xmiss} can be specified as attributes to the input data matrix before invoking attimp (see \code{\link[COSA]{cosa}}).                                             
#'
#'
#' @return 
#' ... STILL NEEDS A DESCRIPTION ... 
#'
#' @include bplot.R
#' @examples
#'   #acquire object at from the attimp function?
#'   #attvalues(X, group[[2]], at$att[1:5])                                         
#'   #attvalues(X, group[[2]], at$att[6:10], all=T)                                  
#'   #attvalues(X, group[[4]], 44, att.names = 'Att 44', group.name = 'Group 4')
#'
#' @seealso
#' \code{\link[COSA]{dendro}}, \code{\link[stats]{hclust}}, \code{\link[COSA]{getclust}}, and \code{\link[COSA]{attvalues}}.
#'
#' @export                                                                               
attvalues <- function(X,group,atts,lX= NULL,att.names=1:ncol(X),         
   group.name='Group', xmiss=NULL,nbins=25,all=F) {                          
   	
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
        	if(class(x) == 'numeric'){
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
   	


   lv=length(atts); lims=matrix(nrow=2,ncol=lv)                                
   if (all) {oldpar=par(mfrow=c(3,lv)); on.exit(par(oldpar))}                 
   else { oldpar=par(mfrow=c(2,lv)); on.exit(par(oldpar))}                     
  
   for (j in 1:lv) { vj=atts[j]                                                
      if (lXt[vj]<4) {                                                         
         z=X[X[,vj]<xmisst,vj]                                                 
         lims[1,j]=min(z); lims[2,j]=max(z)                                    
         m=group[X[group,vj]<xmisst]                                           
         hist(X[m,vj],nclass=nbins,xlim=lims[,j],ylab='Counts',                
            main=paste('Attribute',as.character(att.names[vj])),               
            xlab=group.name)                                                   
      } else {                                                                   
         u=sort(unique(X[,vj]))                                                
         lu=length(u); if (u[lu]>=xmisst) lu=lu-1                              
         bars=rep(0,lu)                                                        
         for (l in 1:lu) { bars[l]=length(group[X[group,vj]==u[l]])}           
         barplot(bars,names.arg=as.character(u[1:lu]),                         
            main=paste('Attribute',as.character(att.names[vj])),               
            xlab=group.name)                                                   
      }                                                                        
   }                                                                           
   for (j in 1:lv) { vj=atts[j]; othr=1:nrow(X); othr=othr[-group]             
      if (lXt[vj]<4) {                                                         
         othr=othr[X[othr,vj]<xmisst]                                          
         hist(X[othr,vj],nclass=nbins,xlim=lims[,j],ylab='Counts',             
            main='',xlab='Others')                                             
      } else {                                                                   
         u=sort(unique(X[,vj]))                                                
         lu=length(u); if (u[lu]>=xmisst) lu=lu-1                              
         bars=rep(0,lu)                                                        
         for (l in 1:lu) { bars[l]=length(othr[X[othr,vj]==u[l]])}             
         barplot(bars,names.arg=as.character(u[1:lu]),                         
            main='',xlab='Others')                                             
      }                                                                        
   }                                                                           
   if (all) {                                                                  
     for (j in 1:lv) { vj=atts[j]                                              
         if (lXt[vj]<4) {                                                      
            z=X[X[,vj]<xmisst,vj]                                              
            hist(z,nclass=nbins,xlim=lims[,j],ylab='Counts',                   
               xlab='All',main='')                                             
         } else {                                                                
            u=sort(unique(X[,vj]))                                             
            lu=length(u); if (u[lu]>=xmisst) lu=lu-1                           
            bars=rep(0,lu)                                                     
            for (l in 1:lu) { bars[l]=length(X[X[,vj]==u[l],vj])}              
            barplot(bars,names.arg=as.character(u[1:lu]),xlab='All',main='')   
         }                                                                     
      }                                                                        
   }                                                                           
   invisible()                                                                 
}
                                                                               
