`opt.RWu` <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for unbiased RW (mu=0)
{
  p0<- mle.rw(y)
  p0<- p0[2]

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.RWu, method="L-BFGS-B", lower=0, control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.RWu, method="BFGS", control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
	cat ("NOTE: ndeps adjusted automatically\n\n")
	cl$ndeps<- p0/1000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.RWu, method="L-BFGS-B", lower=0, control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)
      
    else if (meth=="BFGS")
      w<- try(optim(p0, fn=logL.RWu, method="BFGS", control=cl, hessian=hess, y=y ,pool=pool), silent=TRUE)
	 
	if (class(w)=="try-error")  # if still doesn't work
	  	{ cat ("*")
		  w$par<- NA
		  w$value<- NA }
  }
  
  sv<- -1/diag(w$hessian)
  se<- sqrt(sv)
  w$se<- se
  w$p0<- p0

  # calculate AIC, and AICc (corrected for low n/K)
    names(w$par)<- "vstep"
  K<-1
  w$K<- K
  n<-length(y$mm)-1
  w$AIC<- -2*w$value + 2*K
  w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions
  w$BIC<- -2*w$value + K*log(n)
  
  return (w)
}

