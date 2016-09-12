bmat <- function(dis,dist)
{
	n<-nrow(dist)
	nnn<-matrix(0,nrow=n,ncol=n)
	bb<-matrix(0,nrow=n,ncol=n)
	for (i in 1:n) { nnn[i,i]<-1 }
	for (i in 1:n)
	{
		sb<-0
		for (j in 1:n)
		{
			if (dist[i,j] != 0)
				{
					bb[i,j]<-dis[i,j]/dist[i,j]
					sb<-sb-bb[i,j]
				}
		}
		bb[i,i]<-bb[i,i]+sb
	}
	bmat<- -bb
	return(bmat)
}
