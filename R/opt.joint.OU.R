opt.joint.OU <-
function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize OU model using tree methods
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt) 
 
 ## get initial estimates
 w0<- mle.GRW(x)
 p0<- c(x$mm[1], w0[2]/10, x$mm[length(x$mm)], 0.05)
 names(p0)<- c("anc","vstep","theta","alpha")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, lower=c(NA,0,NA,0.0001), hessian=hess, x=x)
 else 				  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, hessian=hess, x=x) 

 w$p0<- p0
 w$K<- 4
 w$n<- length(x$mm)
 w$AIC <- IC(w, meth="AIC")
 w$AICc<- IC(w, meth="AICc")
 w$BIC <- IC(w, meth="BIC")
 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian))) 
   
 return(w)		
}

