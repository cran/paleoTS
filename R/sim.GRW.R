`sim.GRW` <-
function (ns=20, ms=0, vs=0.1, vp=1, nn=rep(20,ns), tt=1:ns)
# simulates GRW; ns= number of samples, ms=mean and vs=variance of the step distribution,
#  vp=population variance, tt=ages of the samples
{
 MM<- array(dim=ns)
 mm<- array(dim=ns)
 vv<- array(dim=ns)
 dt<- diff(tt)
 
 inc<- rnorm(ns-1, ms*dt, sqrt(vs*dt))	# evolutionary increments
  
 MM<- cumsum(c(0,inc))	# true means
 mm<- MM + rnorm(ns, 0, sqrt(vp/nn))	# true means plus sampling error
 vv<- rep(vp, ns)
 
 gp<- c(ms, vs)
 names(gp)<- c("mstep", "vstep")
 
 res<- as.paleoTS(mm=mm, vv=vv, nn=nn, tt=tt, MM=MM, genpars=gp, label="Created by sim.GRW()") 
 return(res)
}

