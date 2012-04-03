fit.sgs <-
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
    if (!silent) 	cat (i, "\t", GG[,i], "\t", round(w$logL, 3), "\n")
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 if (!silent) cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}
