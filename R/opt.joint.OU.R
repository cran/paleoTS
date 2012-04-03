opt.joint.OU <-
function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize OU model using tree methods
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt) 
 
 ## get initial estimates
 w0<- mle.GRW(x)
 halft<- (x$tt[length(x$tt)]-x$tt[1])/4			# set half life to 1/4 of length of sequence
 p0<- c(x$mm[1], w0[2]/10, x$mm[length(x$mm)], log(2)/halft)
 names(p0)<- c("anc","vstep","theta","alpha")
 #print(p0)
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, lower=c(NA,1e-10,NA,1e-8), hessian=hess, x=x)
 else 				  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, hessian=hess, x=x) 

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='OU', method='Joint', K=4, n=length(x$mm), se=w$se)
  
 return (wc)
}
