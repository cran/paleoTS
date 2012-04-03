opt.Stasis <-
function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.Stasis(y)
  if (p0[2] <= 0 || is.na(p0[2]))	p0[2]<- 1e-7
  names(p0)<- c("theta", "omega")
  if (is.null(cl$ndeps))		cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.Stasis failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}
