`opt.RW.SameVs` <-
function (y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Ms model across multiple sequences
{
  if (class(y)=="paleoTS")
  	stop("Function opt.SameVs() is only meaningful for multiple sequences.\n")
  nseq<- length(y)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	y[[i]]<- pool.var(y[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(y[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(p0m, median(p0v))  # separate Ms, followed by shared Vs for each sequence
  names(p0)<- c(paste("mstep", 1:nseq, sep=""), "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)  
 

  # optimize logL
  ll<- c(rep(NA,nseq), 0)
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, y=y), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, y=y), silent=TRUE)
    else
      w<- try(optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)
	if (class(w)=="try-error")  # if still doesn't work
	  	{   warning("opt.RW.SameVs failed ", immediate.=TRUE)
		    w$par<- NA
		    w$value<- NA }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(-1/diag(w$hessian))
  w$p0<- p0
  w$K<- nseq+1
  n<-0
  for (i in 1:nseq)	n<- n + (length(y[[i]]$mm)-1)
  w$n<- n
  w$AIC<- IC(w, meth="AIC")
  w$AICc<- IC(w, meth="AICc")
  w$BIC<- IC(w, meth="BIC")
  
  return (w)
}

