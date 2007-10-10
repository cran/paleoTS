`mle.GRW` <-
function(y)
# Gives analytical parameter estimates (GRW), assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(y$mm)-1 
 tt<- (y$tt[nn+1]-y$tt[1])/nn
 eps<- 2*pool.var(y)/round(mean(y$nn))  # sampling variance
 dy<- diff(y$mm)
 my<- mean(dy)
 
 mhat<- my/tt
 vhat<- (1/tt)*( (1/nn)*sum(dy^2) - my^2 - eps)
 
 w<- c(mhat, vhat)
 names(w)<- c("mstep", "vstep")
 return(w)
}

