`fit.sgs` <-
function(y, minb=5, oshare=TRUE, pool=TRUE, silent=FALSE, hess=FALSE, meth="L-BFGS-B", model="GRW")
## optimize for stasis-GRW-stasis dynamics (with some min n per section)
{
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)  # pool variances
 ns<- length(y$mm)
 ng<-3
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n\n", "i\tshifts\tlogL\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc) 
  {
    gg<- shift2gg(GG[,i], ns)
	w<- opt.sgs(y, gg, oshare=oshare, hess=hess, meth=meth, model=model)              
    if (!silent) 	cat (i, "\t", GG[,i], "\t", round(w$val, 3), "\n")
    logl[i]<- w$value
    wl[[i]]<- w
  }

 # add more information to results
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ww$n<- ns-1
 ww$shift.start<- GG[,winner]
 ww$AIC<- IC(ww, meth="AIC")
 ww$AICc<- IC(ww, meth="AICc")
 ww$BIC<- IC(ww, meth="BIC")
 ww$all.logl<- logl
 ww$GG<- GG
 if (hess)		ww$se<- sqrt(diag(-1*solve(ww$hessian)))

 if (!silent) cat("\n")
 return(ww)
}

