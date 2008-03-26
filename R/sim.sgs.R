`sim.sgs` <-
function (ns=c(20,20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, times=sum(ns)), tt=1:sum(ns), vp=1)
# simulate stasis-grw-stasis sequence, take theta2 to be final value after grw part
{
  xl<- list()
  for (i in 1:3)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1	
   	 }
   	 
   	 if (i==2)
   	 	xl[[i]]<- sim.GRW(ns[2], ms, vs, nn=nn[start.i:end.i], tt=tt[start.i:end.i], vp=vp)
   	 
   	 else
   	 	xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta, omega=omega, vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }

  ## add offsets
  xl[[2]]$mm<- xl[[2]]$mm + xl[[1]]$MM[ns[1]]
  xl[[3]]$mm<- xl[[3]]$mm + xl[[2]]$MM[ns[2]]
	
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.sgs()"
  y$genpars <- c(theta, omega, ms, vs)
  names(y$genpars)<- c("theta","omega", "ms","vs") 
  return(y)
}

