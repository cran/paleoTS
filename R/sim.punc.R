sim.punc <-
function (ns=c(10,10), theta=c(0,1), omega=rep(0,length(theta)), nn=rep(30,sum(ns)), tt=1:sum(ns), vp=1)
# simulate punctuated sequence; theta and omega are vectors of paramters
# ns is vector of ns in each sub-sequence
{
  nr<- length(theta)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1		
   	 }
   	 
   	 xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta[i], omega=omega[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.punc()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(theta, omega, shft)
  names(y$genpars)<- c(paste("theta",1:nr,sep=""), paste("omega",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
  return (y)  	
}

