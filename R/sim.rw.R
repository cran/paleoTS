"sim.rw" <-
function (ns, sm, sv, vp=1, nn=rep(20,ns), tt=1:ns)
# returns random walk with sampling, based on normal expectation at each time step
#  ns= number of samples, sm=mean and sv=variance of the step distribution,
#  vp=population variance, tt=ages of the samples
{
 MM<- array(dim=ns)
 mm<- array(dim=ns)
 vv<- array(dim=ns)
 dt<- diff(tt)
 
 MM[1]<-0
 x<- rnorm(nn[1], mean=MM[1], sd=sqrt(vp))
 mm[1]<- mean(x)
 vv[1]<- var(x)
 for (i in 2:ns)
 {
    MM[i]<- MM[i-1] + rnorm(1,sm*dt[i-1],sqrt(sv*dt[i-1]))	# prev mean, plus N(ms*t,vs*t)
    x<- rnorm (nn[i], mean=MM[i], sd=sqrt(vp))
    mm[i]<- mean(x)
    vv[i]<- var(x)
 }
 
 res<- as.paleoTS(mm=mm, vv=vv, nn=nn, tt=tt, MM=MM, genpars=c(sm,sv), label="Created by sim.rw()") 
 return(res)
}

