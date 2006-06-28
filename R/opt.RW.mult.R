"opt.RW.mult" <-
function (y, cl=list(fnscale=-1), model=c("RW", "RWu"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates single model across multiple sequences
{
  model<- match.arg(model)
  p0<- mle.rw(y[[1]])

  if (model=="RW")
     ll<- c(NA,0)
  else if (model=="RWu")
     { ll<- 0
       p0<- p0[2]  }

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.mult, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool, model=model), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.mult, method="BFGS", control=cl, hessian=hess, y=y, pool=pool, model=model), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/1000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.mult, method="L-BFGS-B", lower=ll, control=cl, hessian=hess, y=y, pool=pool, model=model), silent=TRUE)
    else if (meth=="BFGS")
      w<- try(optim(p0, fn=logL.mult, method="BFGS", control=cl, hessian=hess, y=y, pool=pool, model=model), silent=TRUE)
	if (class(w)=="try-error")  # if still doesn't work
	  	{   cat ("*")
		    w$par<- NA
		    w$value<- NA }
  }
  
  sv<- -1/diag(w$hessian)
  se<- sqrt(sv)
  w$se<- se
  w$p0<- p0
  if (model=="RW")
     K<- 2
  else if (model=="RWu")
     K<- 1
     
  # calculate AIC, and AICc (corrected for low n/K)
  w$K<- K
  n<-0
  for (j in 1:length(y))
  	n<-n+ length(y[[j]]$mm)-1
  n<- n-length(y)
  w$AIC<- -2*w$value + 2*K
  w$AICc<- w$AIC + (2*K*(K+1))/(n-K-1)  #n is considered to be the number of evolutionary transitions

  return (w)
}

