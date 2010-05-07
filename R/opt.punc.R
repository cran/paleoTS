opt.punc <-
function(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare) 
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
 
  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- 2*ng; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 3*ng-1; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  
  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }
  
  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  names(w$par)<- names(p0)<- pn
  w$p0 <- p0
  w$K <- K
  w$n <- length(y$mm) - 1
  w$AIC<- IC(w, meth="AIC")
  w$AICc<- IC(w, meth="AICc")
  w$BIC<- IC(w, meth="BIC")
  return(w)
}

