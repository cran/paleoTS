`opt.GRW.shift` <-
function(y, ng=2, minb=5, model=1, pool=TRUE, silent=FALSE)
## optimize for shifted GRW dynamics (with some min n per section)
## models:	1  grw (same Vs, diff Ms)
#			2  grw (same Ms, diff Vs)
#			3  urw (diff Vs)
#			4  grw (diff Ms, diff Vs)
{
 ns<- length(y$mm)
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    gg<- shift2gg(GG[,i], ns)
    yl<- split4punc(y, gg)
    #print(yl[[1]]$)
    if (model==1)	{ w<- opt.RW.SameVs(yl, pool=pool);  w$K<- 2*ng; }
    if (model==2)	{ w<- opt.RW.SameMs(yl, pool=pool);  w$K<- 2*ng; }
    if (model>=3)
     {
      wli<- list()
      totS<-0
      totpar<- numeric()
      for (j in 1:ng) 
       {
      	if (model==3){	
      		wli[[j]]<- opt.URW(yl[[j]], pool=pool)
      		names(wli[[j]]$par)<- paste("vstep", j, sep="")   }
      	if (model==4){	
      		wli[[j]]<- opt.GRW (yl[[j]], pool=pool)
      		names(wli[[j]]$par)<- paste(c("mstep","vstep"), j, sep="")   }
      	totS<- totS + wli[[j]]$value
      	totpar<- c(totpar, wli[[j]]$par)
       }
      w<- list(value=totS, par=totpar, K=(model-1)*ng -1)	
     }
         
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
 # hessian not easily recoverable in this optimization
 #if (hess)		ww$se<- sqrt(diag(-1*solve(ww$hessian)))

 if (!silent) cat("\n")
 return(ww)
}

