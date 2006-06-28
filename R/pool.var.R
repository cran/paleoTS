`pool.var` <-
function(y,nn=NULL, ret.paleoTS=FALSE)
# calc pooled variance
# y is either a paleoTS object, or an array of variances
{
 # check if y is a paleoTS object
 if (class(y)=="paleoTS")
    vp<- sum(y$vv*(y$nn-1))/sum(y$nn-1)
 else
    vp<- sum(y*(nn-1))/sum(nn-1)
    
 if (ret.paleoTS)
 {
   yn<- y
   yn$vv<- rep(vp, length(y$mm))
   return(yn)		
 }
 else { return (vp)}
}

