sim.GRW.shift <-
function (ns=c(10,10), ms=c(0,1), vs=c(0.5,0.5), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate sequence of GRW with different parameter values in different segments
# ns is vector of ns in each sub-sequence, similar for other parameters
{
  nr<- length(ms)
  cns<- cumsum(ns)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- cns[i-1]+1
   	  end.i<- cns[i-1]+ns[i]	
   	 }
   	 
   	 xl[[i]]<- sim.GRW(ns=ns[i], ms=ms[i], vs=vs[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   	 if (i>1)
   	 	xl[[i]]$mm<- xl[[i]]$mm + xl[[i-1]]$mm[length(xl[[i-1]]$mm)]
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.GRW.shift()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(ms, vs, shft)
  names(y$genpars)<- c(paste("mstep",1:nr,sep=""), paste("vstep",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))  
  return (y)  	
}
