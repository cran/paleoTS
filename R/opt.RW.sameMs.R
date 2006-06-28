"opt.RW.sameMs" <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Ms model across multiple sequences
{
  # generate initial parameter estimates {ms, v1,..vk}
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
  p0<- c(median(p0m), p0v)  # shared Ms, followed by separate Vs for each sequence
  	 

  # optimize logL
  ll<- c(NA,0)
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.sameMs, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.sameMs, method="BFGS", control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/1000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.sameMs, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
    else if (meth=="BFGS")
      w<- try(optim(p0, fn=logL.sameMs, method="BFGS", control=cl, hessian=hess, y=y, pool=pool), silent=TRUE)
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
     
  # calculate AIC, and AICc (corrected for low n/K)
  n<-0
  for (j in 1:length(y))
  	n<-n+ length(y[[j]]$mm)-1
  #n<- n-length(y)   # this line shouldn't be there!
  w$AIC<- -2*w$value + 2*K
  w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions

  return (w)
}

