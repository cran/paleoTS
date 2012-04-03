opt.URW <-
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
	cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
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
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='AD', K=1, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}
