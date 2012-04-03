opt.RW.SameMs <-
function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Ms model across multiple sequences
{
  if (class(yl)=="paleoTS")
  	stop("Function opt.SameMs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(yl[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(median(p0m), p0v)  # shared Ms, followed by separate Vs for each sequence
  names(p0)<- c("mstep", paste("vstep", 1:nseq, sep=""))
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)

  # optimize logL
  ll<- c(NA, rep(0,nseq))
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
    else 
      w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)
	#if (class(w)=="try-error")  # if still doesn't work
	#  	{   warning("opt.RW.SameMs failed ", immediate.=TRUE)
	#	    w$par<- NA
	#	    w$value<- NA }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameMs.Mult", method='AD', K=nseq+1, n=n, se=w$se)
  
  return (wc)
}
