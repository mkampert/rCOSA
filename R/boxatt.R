#' Barplots of univariate attribute dispersion within groups
#'
#' This function plots the attribute weights over groups.
#' 
#' @param data data matrix obtained from COSA analysis
#' @param imps typically a vector with attribute indices ordered on their importance
#' @param grps ...
#' @param main ...
#' @param ylab ...
#' @param xlab ...
#' @param SCALE logical TRUE/FALSE, whether to scale the data
#' @param colors  ... 
#' @param ... named arguments to be passed to the default method.
#' 
#'
#' @return 
#' values of attribute importances in requested order.                         
#'
#' @examples
#'  make sure an at object is available...
#'
#' @seealso
#' \code{\link[COSA]{attimp}}
#' 
#' 
#' @export                                                                               


boxatt <- function(data, imps, grps, 
                   main = "Consider a maintitle here", 
                   ylab = "Consider a label here", xlab = 'Consider a label here',
                   SCALE = TRUE, colors = palette(), ...){
  
  if(SCALE) p <- autoscale(data)
  boxplot(p[grps[[1]] ,imps],
          main = main, xlab = xlab, ylab = ylab,
          outline = FALSE, border = colors[1], col = colors[1],
          names = 1:length(imps), ...)
  
  colMdns <- matrix(NA, ncol = length(imps), nrow = length(grps))
  for(i in 1:length(grps)){
    colMdns[i,] <- apply(p[grps[[i]], imps], 2, median, na.rm = TRUE, ...)
    text(colMdns[i,], col= colors[i + 1], label='-')
  }
}



# Examples
#boxatt(data = ApoE3, imps = attWild$att[1:100], 
#       grps =  list(1:18, 19:38), 
#       main = "Wildtype (red) vs Transgenic (blue)", 
#       ylab = 'Scaled attribute values', 
#       xlab = '100 most important attributes',
#       colors = c('red','white','blue'))