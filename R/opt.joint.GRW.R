opt.joint.GRW <-
function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize GRW model using alternate formulation
{
 ## check if pooled, make start at tt=0
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt)
  
 ## get initial estimates
 p0<- array(dim=3)
 p0[1]<- x$mm[1]	
 p0[2:3]<- mle.GRW(x)
 if (p0[3]<=0)	p0[3]<- 1e-7
 names(p0)<- c("anc", "mstep", "vstep")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
 if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, lower=c(NA,NA,0), hessian=hess, x=x)
 else 					w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, hessian=hess, x=x)

 w$p0<- p0
 w$K<- 3
 w$n<- length(x$mm)
 w$AIC <- IC(w, meth="AIC")
 w$AICc<- IC(w, meth="AICc")
 w$BIC <- IC(w, meth="BIC")
 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian))) 

 return(w)		
}

