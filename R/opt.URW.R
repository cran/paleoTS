`opt.URW` <-
function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.URW(y)
  if (p0 <= 0)	p0<- 1e-7
  names(p0)<- "vstep"
  if (is.null(cl$ndeps))		cl$ndeps<- p0/1e4
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.URW failed ", immediate.=TRUE)
		  w$par<- NA
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(-1/diag(w$hessian))
  w$p0<- p0
  w$K<- 1
  w$n<- length(y$mm)-1
  w$AIC<- IC(w, meth="AIC")
  w$AICc<- IC(w, meth="AICc")
  w$BIC<- IC(w, meth="BIC")
  
  return (w)
}

