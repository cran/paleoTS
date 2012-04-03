opt.covTrack.Mult <-
function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
# y and z are lists of paleoTS, and covariates, respectively
{
 if (class(yl) == "paleoTS") 
      { stop("opt.track.Mult is only for multiple paleoTS sequences\n") }
 nseq <- length(yl)
 if (pool) {
        for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE) }
        
 # check lengths of z
 for (i in 1:nseq)
 {
    ns<- length(yl[[i]]$mm)
    if(length(zl[[i]])==length(yl[[i]]$mm))	
  	 { 
  	  zl[[i]]<- diff(zl[[i]])
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(zl[[i]]) != ns-1)  stop("Covariate length [", length(zl[[i]]), "] does not match the sequence length [", ns, "]\n" )
 	
 }

        
 p0<- c(0,var(yl[[1]]$mm))  # lousy method for initial guess!
 names(p0)<- c("b", "evar")
 w<- optim(p0, fn=logL.Mult.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,0), hessian=hess, yl=yl, zl=zl)
 
 if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
 else		w$se<- NULL
 
 ff<- function(x) length(x$mm)-1
 n<- sum(sapply(yl, ff))
 
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='trackCovariate.Mult', method='AD', K=2, n=n, se=w$se)
 return(wc)	
}
