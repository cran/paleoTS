sim.covTrack <-
function (ns=20, b=1, evar=0.1, z, nn=rep(20, times=ns), tt=0:(ns-1), vp=1)
# simulates tracking optimum, with covariate z
{
  # check on length of covariate
  if(length(z)==ns)	
  	{ 
  	 z<- diff(z)
	 warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
  
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- rnorm(ns-1, b*z, sqrt(evar))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))  # add sampling error
  vv <- rep(vp, ns)
  gp <- c(b, evar)
  names(gp) <- c("slope.b", "evar")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
        genpars = gp, label = "Created by sim.covTrack()")
  return(res)	
}
