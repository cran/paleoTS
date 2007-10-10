`opt.GRW` <-
function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.GRW(y)
  if (p0[2] <= 0)	p0[2]<- 1e-7
  names(p0)<- c("mstep", "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.GRW failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(-1/diag(w$hessian))
  w$p0<- p0
  w$K<- 2
  w$n<- length(y$mm)-1
  w$AIC<- IC(w, meth="AIC")
  w$AICc<- IC(w, meth="AICc")
  w$BIC<- IC(w, meth="BIC")
  
  return (w)
}

