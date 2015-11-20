autoscale=function(x, xmi=NA, xmo=NA, robust=F) 
{                              
  y<-matrix(nrow=nrow(x),ncol=ncol(x))                                         
  if (!is.na(xmi)) x[x==xmi]<-NA                                               
  for (i in 1:ncol(x)) {                                                      
    u<-x[!is.na(x[,i]),i]                                                     
    if (robust) { scl<-IQR(u)/1.349;                                          
                  if (scl > 0) u<-(u-median(u))/scl
                  if (scl == 0) u<-(u-mean(u))/sqrt(var(u))                                      
    }                                                                        
    else { scl<-var(u);                                                       
           if (scl > 0) { scl<-sqrt(scl); u<-(u-mean(u))/scl}                      
    }                                                                        
    y[!is.na(x[,i]),i]<-u; y[is.na(x[,i]),i]<-xmo                              
  }                                                                           
  invisible(y)                                                                
}