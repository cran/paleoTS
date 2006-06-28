`opt.RW` <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for general RW (mu is estimated)
# y is an paleoTS object
{
  p0<- mle.rw(y)
  if (p0[2] <= 0)	p0[2]<- 1e-5

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.RW, method="L-BFGS-B", lower=c(NA,0), control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)
  else if (meth=="BFGS")
   w<- try(optim(p0, fn=logL.RW, method="BFGS", control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	  cl$ndeps<- p0/1000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.RW, method="L-BFGS-B", lower=c(NA,0), control=cl, hessian=hess, y=y,pool=pool), silent=TRUE)
    else if (meth=="BFGS")
     w<- try(optim(p0, fn=logL.RW, method="BFGS", control=cl, hessian=hess, y=y ,pool=pool), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  cat ("@")
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  sv<- -1/diag(w$hessian)
  se<- sqrt(sv)
  if (hess)
  	w$se<- se
  else w$se<- NULL
  w$p0<- p0

  # calculate AIC, and AICc (corrected for low n/K)
  names(w$par)<- c("mstep", "vstep")
  K<-2
  w$K<- K
  n<-length(y$mm)-1
  w$AIC<- -2*w$value + 2*K
  w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions
  w$BIC<- -2*w$value + K*log(n)
  
  return (w)
}

