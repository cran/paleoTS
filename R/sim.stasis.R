`sim.Stasis` <-
function(ns=20, theta=0, omega=0, vp=1, nn=rep(20,ns), tt=1:ns)
# simulate stasis
{
 xmu<- rnorm(ns, mean=theta, sd=sqrt(omega))
 xobs<- xmu + rnorm(ns, 0, sqrt(vp/nn))
 gp<- c(theta, omega)
 names(gp)<- c("theta", "omega")
   	
 x <- as.paleoTS(mm=xobs,vv=rep(vp,ns),nn=nn,tt=tt,MM=xmu,genpars=gp,label="Created by sim.stasis()") 
 return(x)	 	 	
}

