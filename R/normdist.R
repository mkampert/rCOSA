#' Normalize the distance 
#'
#'
#' Normalize the dissimilarities sucht that the squared norm of the Dissimilarities equals N, the number of objects.
#'
#'
#' @param D A matrix of distances of class dist, or a vector or matrix containing distances.
#'
#'
#' @return
#' Returns a normalised version of D of similar class and mode.


normdist=function(D)
{
        fac<-sum(D**2)
        n<-dimdata(D)
        D<-D*((n/fac)**.5)
        return(D)
}
