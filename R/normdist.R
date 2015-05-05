normdist=function(dd)
{
        fac<-sum(dd**2)
        n<-dimdata(dd)
        dd<-dd*((n/fac)**.5)
        return(dd)
}
