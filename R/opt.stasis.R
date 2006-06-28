`opt.stasis` <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# optimize estimated of stasis model
{
 p0<- c(mean(y$mm), var(y$mm))
 K<-2
 cl$ndeps<- p0/1000
 
 if(meth=="L-BFGS-B")	w<- optim(p0, fn=logL.stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y, pool=pool)
 else			 		w<- optim(p0, fn=logL.stasis, method=meth, control=cl, hessian=hess, y=y, pool=pool)

 sv<- -1/diag(w$hessian)
 se<- sqrt(sv)
 if (hess)
  	w$se<- se
 else w$se<- NULL
 w$p0<- p0

 # calculate AIC, and AICc (corrected for low n/K)
 names(w$par)<- c("theta", "omega")
 w$K<- K
 n<-length(y$mm)
 w$AIC<- -2*w$value + 2*K
 w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions
 w$BIC<- -2*w$value + K*log(n)
  
 return (w) 	
}

