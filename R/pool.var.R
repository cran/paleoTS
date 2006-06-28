"pool.var" <-
function(y,nn=NULL)
# calc pooled variance
# y is either a paleoTS object, or an array of variances
{
 # check if y is a paleoTS object
 if (class(y)=="paleoTS")
    vp<- sum(y$vv*(y$nn-1))/sum(y$nn-1)
 else
    vp<- sum(y*(nn-1))/sum(nn-1)
    
 return (vp)
}

