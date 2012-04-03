opt.covTrack <-
function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE) 
{
    # check if z is of proper length; first difference if necessary
    ns<- length(y$mm)
    if(length(z)==length(y$mm))	
  	 { 
  	  z<- diff(z)
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(z) != ns-1)  stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n" )

    
    # get initial estimates by regression
    reg<- lm(diff(y$mm) ~ z-1)
    p0<- c(coef(reg), var(resid(reg)))
    names(p0) <- c("b", "evar")
    
    # pool variances if needed and do optimization
    if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
    if (is.null(cl$ndeps)) 
        cl$ndeps <- abs(p0/10000)
    if (meth == "L-BFGS-B") 
        w <- optim(p0, fn=logL.covTrack, method = meth, lower = c(NA, 0), control = cl, hessian = hess, y=y, z=z)
    else  w<- optim(p0, fn=logL.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)
   

    if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
    else w$se <- NULL
    
	wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='AD', K=2, n=length(y$mm)-1, se=w$se)
}
