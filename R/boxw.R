#' Barplots of attribute weights within groups
#'
#' This function plots the attribute weights over groups.
#' 
#' @param W weights matrix obtained from COSA analysis
#' @param grpnr object row numbers of a specific group 
#' @param attr attribute column numbers ordererd from important to less important. To find out the order use \code{attimp}.   
#' @param ... named arguments to be passed to the default method.
#' 
#'
#' @return 
#' values of attribute importances in requested order.                         
#'
#' @examples
#'  #make sure an at object is available...
#'
#' @seealso
#' \code{\link[COSA]{attimp}}
#' 
#' 
#' @export                                                                               


boxw <- function(W, grpnr, attr, ...){
  boxplot( W[grpnr, attr], ...)
}

