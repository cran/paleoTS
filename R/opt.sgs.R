opt.sgs <-
function(y,gg,cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare=TRUE, model="GRW")
# do optimization for sgs model
{
 yl<- split4punc(y, gg)
 hat<- mle.GRW(yl[[2]])
 if (hat[1] == 0)	hat[1]<- 1e-3
 if (hat[2] < 1e-4)	hat[2]<- 1e-4
 if (model=="URW")	p0<- c(hat[2], mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 else 				p0<- c(hat, mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 
 if (oshare)	
 	{ p0<- append(p0, mean(var(yl[[1]]$mm), var(yl[[3]]$mm)))
 	  if (model=="GRW")		{ K<-7; pn<- c("ms","vs","theta1","theta2","omega"); lw<- c(NA,0,NA,NA,0) }
 	  else					{ K<-6; pn<- c("vs","theta1","theta2","omega"); lw<- c(0,NA,NA,0) }	}  
 else
 	{ p0<- append(p0, c(var(yl[[1]]$mm), var(yl[[3]]$mm)) )
 	  if (model=="GRW")	{ K<- 8; pn<- c("ms","vs","theta1","theta2","omega1","omega2"); lw<- c(NA,0,NA,NA,0,0) }
 	  else 				{ K<- 7; pn<- c("vs","theta1","theta2","omega1","omega2"); lw<- c(0,NA,NA,0,0) }
 	}
 
 cl$ndeps <- p0/100
 #cat (p0, "\n")
 if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.sgs.omega, gg=gg, method=meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn=logL.sgs.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.sgs, gg=gg, method= meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn = logL.sgs, gg=gg, method = meth, control = cl, hessian = hess, y = y)
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

