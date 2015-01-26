#' Title for this function
#'
#'    Description of the function here                                         
#'
#' @param x 
#' @param xlab...

bplot <- function(x,xlab='',ylab='',zlim=c(0,1),names=NULL,horiz=F,trans=1,ntk=20,
   names.size=1){                                                             

   if (trans==1) {                                                             
      if (horiz) {                                                             
         barplot(x,names=names,horiz=T,xlab=xlab,ylab=ylab,xlim=zlim,las=1,    
            cex.names=names.size)                                              
      }                                                                        
      else {                                                                   
         barplot(x,names=names,horiz=F,xlab=xlab,ylab=ylab,ylim=zlim,las=2,    
            cex.names=names.size)                                              
      }                                                                        
   }                                                                           
   else if (horiz) {                                                           
      barplot(x**trans,names=names,horiz=T,xlab=xlab,ylab=ylab,las=1,          
         xlim=zlim**trans,axes=F,cex.names=names.size)                         
      lab <- pretty(x,n=ntk)                                                   
      axis(side = 1, at = lab**trans, labels = as.character(lab))              
   }                                                                           
   else {                                                                      
      barplot(x**trans,names=names,xlab=xlab,ylab=ylab,ylim=zlim**trans,axes=F,
        las=2,cex.names=names.size)                                            
      lab <- pretty(x,n=ntk)                                                   
      axis(side = 2, at = lab**trans, labels = as.character(lab))              
   }                                                                           
}
