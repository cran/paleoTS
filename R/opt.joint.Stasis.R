opt.joint.Stasis <-
function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize Stasis model using alternate formulation
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt)
  
 ## get initial estimates
 p0<- mle.Stasis(x)
 if(p0[2]<=0 || is.na(p0[2]))	p0[2]<- 1e-7
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-9
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, lower=c(NA,0), hessian=hess, x=x)
 else				  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, hessian=hess, x=x)
 
 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='Joint', K=2, n=length(x$mm), se=w$se)
  
 return (wc)
	
}
