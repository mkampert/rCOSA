#' Hierarchical Clustering
#'
#' Performs a hierarchical cluster analysis on a set of dissimilarities using the R procedure \code{\link[stats]{hclust}} and \code{\link[stats]{dendrogram}}. 
#' 
#' 
#' @param D dissimilarities/distances of class dist.
#' @param main plot title
#' @param leaflab labels of the leaves, default is set to 'none'
#' @param cex.main redundant argument
#' @param denplot logical, whether to plot dendrogram. Default is set to \code{TRUE}.
#' @param sub subtitle
#' @param xlab label horizontal axis
#' @param ylab label vertical axis
#' @param the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of \code{"ward"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#' @param groupnr when known, each object can be given a number to which group it belongs
#' @param colv vector containing color names for each group number
#' @param ... arguments which can be parsed to either the \code{\link[mva]{hclust}} and \code{\link[stats]{dendrogram}} procedures.
#' @return 
#' Returns an invisible object of class dendrogram. If denplot is set to \code{TRUE} then the dendrogram is also plotted, and  when the groups are known one can ask for a colored dendrogram.
#' @export                                                                               
hierclust <- function(D, main = "Hierarchical Clustering", leaflab = c("perpendicular", "textlike", "none")[3], cex.main = NULL, denplot = TRUE, sub = "", xlab = "", ylab = "", method = 'average', groupnr = NULL, colv = palette()[c(8,2,4,3,5,6,7,1)], ...){

	# This is code based on the hclust function
	hc <- hclust(D, method = method, ...)
	hc$dendro <- as.dendrogram(hc, ...)

	if(denplot){
		plot(hc$dendro,sub = sub, xlab = xlab, ylab = ylab, leaflab = leaflab, main = main, ...)
		
		if(length(groupnr) == attributes(hc$dendro)$members){
			cols <- colv[groupnr + 1]
			for(i in 1:length(groupnr)){
				indx <- (hc$merge[,1] == -hc$ord[i]) | (hc$merge[,2] == -hc$ord[i])
				segments(i, 0,i, hc$height[indx], col = cols[hc$ord[i]])
			}
		}
		
	}
	invisible(hc)
}                                     
                                                                               
