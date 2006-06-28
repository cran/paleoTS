`opt.RW.sameVs` <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Vs model across multiple sequences
{
  # generate initial parameter estimates {m1,..mk, vs}
  nseq<- length(y)
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.rw(y[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-5	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(p0m, median(p0v))  # separate Ms, followed by shared Vs
  	 

  # optimize logL
  ll<- c(rep(NA,nseq),0)
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.sameVs, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.sameVs, method="BFGS", control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/1000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.sameVs, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
    else if (meth=="BFGS")
      w<- try(optim(p0, fn=logL.sameVs, method="BFGS", control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
	if (class(w)=="try-error")  # if still doesn't work
	  	{   cat ("*")
		    w$par<- NA
		    w$value<- NA }
  }
  
  sv<- -1/diag(w$hessian)
  se<- sqrt(sv)
  w$se<- se
  w$p0<- p0
  K<- nseq+1
  names(w$par)<- c(paste("mstep", 1:nseq, sep=""), "vstep")
     
  # calculate AIC, and AICc (corrected for low n/K)
  n<-0
  for (j in 1:length(y))
  	n<-n+ length(y[[j]]$mm)-1
  #n<- n-length(y)   # this line shouldn't be there!
  w$AIC<- -2*w$value + 2*K
  w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions
  w$BIC<- -2*w$value + K*log(n)

  return (w)
}

