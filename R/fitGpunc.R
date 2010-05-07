fitGpunc <-
function(y, ng=2, minb=5, pool=TRUE, oshare=TRUE, silent=FALSE, hess=FALSE, ...)
## optimize punctuation models (with some min n per section)
{
 if(ng==1)  # if only one grouping, same as stasis model
 {
   warning('Fitting stasis model (because ng=1)')
   ww<- opt.Stasis(y, pool=pool, hess=hess, ...)
   return(ww)	
 }
 
 ns<- length(y$mm)
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("Total # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    gg<- shift2gg(GG[,i], ns)
    w<- opt.punc(y, gg, oshare=oshare, pool=pool, hess=hess, ...)
    logl[i]<- w$value
    wl[[i]]<- w
  }
 cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ww$shift.start<- GG[,winner]
 ww$all.logl<- logl
 ww$GG<- GG
 if (oshare) pn<- c(paste("theta",1:ng,sep=""), "omega")
 else 		 pn<- c(paste("theta",1:ng,sep=""), paste("omega",1:ng,sep=""))
 names(ww$par)<- pn
 return(ww)
}

