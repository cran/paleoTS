"mle.rw" <-
function(x)
# Gives analytical parameter estimates, assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(x$mm)-1 
 tt<- (x$tt[nn+1]-x$tt[1])/nn
 eps<- 2*pool.var(x)/round(mean(x$nn))  # sampling variance
 dx<- diff(x$mm)
 mx<- mean(dx)
 
 mhat<- mx/tt
 vhat<- (1/tt)*( (1/nn)*sum(dx^2) - mx^2 - eps)
 
 w<- c(mhat, vhat)
 names(w)<- c("ms.hat", "vs.hat") 
 return(w)
}

