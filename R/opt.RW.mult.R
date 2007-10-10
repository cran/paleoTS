`opt.RW.Mult` <-
function (y, cl=list(fnscale=-1), model=c("GRW", "URW"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates single model across multiple sequences
# pool=TRUE will pool variances _within_ sequences
{
  if (class(y)=="paleoTS")
     stop("opt.RW.mult is onlt for multiple paleoTS sequences\n")
  nseq<- length(y)
     
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	y[[i]]<- pool.var(y[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
     	ll<- c(NA,0)
     	p0<- mle.GRW(y[[1]])	
     	K<- 2	 	}
  else if (model=="URW")
     { ll<- 0
       p0<- mle.URW(y[[1]])
       K<- 1	  	}
  
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, y=y, model=model), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, y=y, model=model), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, y=y, model=model), silent=TRUE)
    else 
      w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, y=y, model=model), silent=TRUE)
	if (class(w)=="try-error")  # if still doesn't work
	  	{   warning("opt.RW.Mult failed ", immediate.=TRUE)
		    w$par<- NA
		    w$value<- NA }
  }
  
  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(-1/diag(w$hessian))
  w$p0<- p0
  w$K<- K
  n<-0
  for (i in 1:nseq)	n<- n + (length(y[[i]]$mm)-1)
  w$n<- n
  w$AIC<- IC(w, meth="AIC")
  w$AICc<- IC(w, meth="AICc")
  w$BIC<- IC(w, meth="BIC")

  return (w)
}

