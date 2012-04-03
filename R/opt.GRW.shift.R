opt.GRW.shift <-
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
    #print(yl[[1]])
    if (model==1)	{ w<- opt.RW.SameVs(yl, pool=pool); w$modelName<- paste('GRWsameVs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model==2)	{ w<- opt.RW.SameMs(yl, pool=pool); w$modelName<- paste('GRWsameMs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model>=3)
     {
      wli<- list()
      totS<-0
      totpar<- numeric()
      for (j in 1:ng) 
       {
      	if (model==3){	
      		wli[[j]]<- opt.URW(yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste("vstep", j, sep="")   }
      	if (model==4){	
      		wli[[j]]<- opt.GRW (yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste(c("mstep","vstep"), j, sep="")   }
      	totS<- totS + wli[[j]]$logL
      	totpar<- c(totpar, wli[[j]]$parameters)
      	kk<- length(totpar)+ng-1
       }
    
      ifelse(model==3, mn<- 'URW-shift', mn<- 'GRW-shift')
      w<- as.paleoTSfit(logL=totS, parameters=totpar, modelName=paste(mn, ng-1, sep='-'), method='AD', K=kk, n=ns-1, se=NULL)	
     }
         
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
 
}
