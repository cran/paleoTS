fitGpunc <-
function(y, ng=2, minb=5, pool=TRUE, oshare=TRUE, method=c('AD', 'Joint'), silent=FALSE, hess=FALSE, ...)
## optimize punctuation models (with some min n per section)
{
 method<- match.arg(method) 
 
 if(ng==1)  # if only one grouping, same as stasis model
 {
   warning('Fitting stasis model (because ng=1)')
   if (method=='AD')	ww<- opt.Stasis(y, pool=pool, hess=hess, ...)
   else if (method=='Joint')	ww<- opt.joint.Stasis(y, pool=pool, hess=hess, ...)
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
    # the different gg for AD and J is required for the interpretation of the "shift" parameters to be the same across parameterizations
    ggA<- shift2gg(GG[,i], ns)
    ggJ<- shift2gg(GG[,i]+1, ns)
    if(method=='AD')    		 w<- opt.punc(y, ggA, oshare=oshare, pool=pool, hess=hess, ...)
    else if (method=='Joint')    w<- opt.joint.punc(y, ggJ, oshare=oshare, pool=pool, hess=hess, ...)
    logl[i]<- w$logL
    wl[[i]]<- w
  }
 cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}
