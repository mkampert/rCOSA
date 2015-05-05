rotat <- function(xx)
{     rot<-eigen((t(xx))%*%xx)
      rot<-rot[[2]]
      xx<-xx%*%rot
      return(xx)
}
