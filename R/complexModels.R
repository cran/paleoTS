##### Punctuations functions #####


# internal function, not exported
cat.paleoTS<- function (y)
# concatenates multiple paleoTS objects, with y a list of paleoTS objects
{
 x<- y[[1]]
 for (i in 2:length(y))
   {
   	x$mm<- append(x$mm, y[[i]]$mm)
   	x$vv<- append(x$vv, y[[i]]$vv)
   	x$tt<- append(x$tt, y[[i]]$tt)
   	x$MM<- append(x$MM, y[[i]]$MM)
   	x$nn<- append(x$nn, y[[i]]$nn)
   }

  return (x)
}


#' Simulate a punctuated time-series
#' @description Simulates punctuated trait evolution with punctuations that are rapid relative
#' to the spacing of samples. In practice, the time-series is divided into two or more
#' segments, each of which has its own mean and variance.
#'
#' @param ns vector of the number of samples in each segment
#' @param theta vector of means, one for each segment
#' @param omega vector of variances, one for each segment.
#' @param nn vector of sample sizes, one for each population
#' @param tt vector of times (ages), one for each population
#' @param vp phenotypic variance within each population
#'
#' @details Segments are separated by punctuations. Population means in the ith segment are
#' drawn randomly from a normal distribution with a mean equal to ith element of \code{theta}
#' and variance equal to the ith element of \code{omega}. The magnitudes of punctuations are
#' determined by the differences in adjacent \code{theta} values.
#'
#' @return a \code{paleoTS} object with the simulated time-series.
#' @seealso \code{\link{fitGpunc}}
#' @export
#'
#' @examples
#' x <- sim.punc(ns = c(15, 15), theta = c(0,3), omega = c(0.1, 0.1))
#' plot(x)
sim.punc<- function (ns=c(10,10), theta=c(0,1), omega=rep(0,length(theta)), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
{
  nr<- length(theta)
  xl<- list()
  for (i in 1:nr)
   {
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1
   	 }

   	 xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta[i], omega=omega[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }

  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.punc()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(theta, omega, shft)
  names(y$genpars)<- c(paste("theta",1:nr,sep=""), paste("omega",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
  return (y)
}

# internal function, not exported
split4punc<- function (y, gg, overlap=TRUE)
# divides a paleoTS object (y) into a several paleoTS objects, according to vector 'gg'
# gg is a vectors of 1,2,3.. indicating groupings
# overlap=TRUE means that the adjacent samples are included
{
  yl<- list()
  ng<- max(gg)
  for (i in 1:ng)
   {
   	 ok<- gg==i
   	 if(i>1 & overlap==TRUE)	ok[max(which(gg==i-1))]<- TRUE   # this is now right!

   	 yl[[i]]<- as.paleoTS(y$mm[ok],y$vv[ok],y$nn[ok],y$tt[ok],y$MM[ok])
   }
  return (yl)
}


# internal function, not exported
logL.punc<- function (p, y, gg)
# logL of punctuation, with shifts
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]

  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S+ logL.Stasis(p=c(th[i], om[i]), xl[[i]])

  return (S)
}

# internal function, not exported
logL.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]

  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S + logL.Stasis(p=c(th[i], om), xl[[i]])

  return (S)
}


#' Fit a model of trait evolution with specified punctuation(s)
#'
#' @param y a \code{paleoTS} object
#' @param gg vector of indices indicating different segments
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param meth optimization algorithm, passed to \code{optim()}
#' @param hess if TRUE, return standard errors of parameter estimates from the
#' @param oshare logical, if TRUE, variance assumed to be shared (equal) across segments
#'
#' @details The sequence is divided into segments, which are separated by punctuations. Means for
#' each segment are given by the vector \code{theta} with variances given by the vector
#' \code{omega} (or a single value if \code{oshare = TRUE}). This function calls \code{optim} to numerically fit this model to a time-series, y.
#'
#' @note These functions would be used in the uncommon situation in which there
#'   is a prior hypothesis as to where the punctuation(s) take place.  Normally
#'   users will instead use the function \code{fitGpunc}, which uses these
#'   functions to fit a range of possible timings for the punctuations.
#'
#' @seealso \code{\link{fitGpunc}}
#'
#' @return a \code{paleoTSfit} object with the results of the model fitting
#' @export
#'
#' @examples
#' x <- sim.punc(ns = c(15, 15), theta = c(0,3), omega = c(0.1, 0.1))
#' w.sta  <- fitSimple(x, model = "Stasis", method = "Joint")
#' w.punc <- opt.joint.punc(x, gg = rep(1:2, each = 15), oshare = TRUE)
#' compareModels(w.sta, w.punc)
opt.punc<- function(y, gg, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare)
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences

  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- ng + 1; pn<- c(pn, "omega")}
  else				  { p0 <- c(mg, mv); K <- 2*ng; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn

  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }

  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)

  return(wc)
}

#internal function, not exported
logL.joint.punc<- function (p, y, gg)
# logL of punctuation, with shifts
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]

  M<- th[gg]  # vector of MVN means
  VV<- diag(om[gg] + y$vv/y$nn)  # vcv matrix
  #S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  S<- mnormt::dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

  return (S)
}

# internal function, not exported
logL.joint.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]

  M<- th[gg]  # vector of MVN means
  omv<- rep(om, max(gg))
  VV<- diag(omv[gg] + y$vv/y$nn)  # vcv matrix
  #S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  S<- mnormt::dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

  return (S)
}

#' @describeIn opt.punc  fits the punctuation model using the joint parameterization
#' @export
opt.joint.punc<- function(y, gg, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare)
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- ng + 1; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 2*ng; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn

  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.joint.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.joint.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }

  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)
}


# internal function, not exported
shifts<- function(ns, ng, minb=7)
{
	aaL<- utils::combn(ns, ng-1, simplify=FALSE)
	aa<- matrix(unlist(aaL), nrow=ng-1)
	top<- rep(0, ncol(aa))
	bot<- rep(ns, ncol(aa))
	aaM<- unname(rbind(top, aa, bot))

	daa<- apply(aaM,2,diff)
	ok<- apply(daa, 2, function(x) all(x >= minb))

	if(ng>2)	return(aa[,ok])
	else		return(matrix(aa[,ok], nrow=1))
}


# internal function, not exported
shift2gg<- function (ss, ns)
# ss is vector of shift points, ns is # samples
{
  z<- c(0,ss,ns+1)
  cc<- cut(1:ns, breaks=z, right=TRUE)
  gg<- as.numeric(cc)
  return(gg)
}



#' Fit trait evolution model with punctuations estimated from the data
#'
#' @param y a \code{paleoTS} object
#' @param ng number of groups (segments) in the sequence
#' @param minb minimum number of populations within each segment
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param oshare logical, if TRUE, variance assumed to be shared (equal) across segments
#' @param method parameterization to use: \code{Joint} or \code{AD}; see Details
#' @param silent logical, if TRUE, progress updates are suppressed
#' @param hess if TRUE, standard errors computed from the Hessian matrix are returned
#' @param parallel logical, if TRUE, the analysis is done in parallel
#' @param ... other arguments, passed to optimization functions
#'
#' @details This function tests all possible shift points for punctuations, subject to the
#' constraint that the number of populations in each segment is always >= \code{minb}. The
#' shiftpoint yielding the highest log-likelihood is returned as the solution, along with
#' the log-likelihoods (\code{all.logl}) of all tested shift points (\code{GG}). \cr \cr
#'
#' The function uses \code{opt.punc} (if \code{method = "AD"}) or \code{opt.joint.punc}
#' (if \code{method = "Joint"}) to do the fitting.
#'
#' @note
#' Calculations can be speeded up by setting \code{parallel = TRUE}, which uses functions from
#' the \code{\link{doParallel}} package to run the bootstrap replicates in parallel, using
#' one fewer than the number of detected cores.
#'
#' @return  a \code{paleoTSfit} object with the results of the model-fitting.
#' @seealso \code{\link{fit9models}}, \code{\link{sim.punc}}
#' @export
#' @import foreach
#' @import doParallel
#'
#' @examples
#' x <- sim.punc(ns = c(15, 15), theta = c(0,3), omega = c(0.1, 0.1))
#' w.punc <- fitGpunc(x, oshare = TRUE)
#' plot(x, modelFit = w.punc)
fitGpunc<-function (y, ng = 2, minb = 7, pool = TRUE, oshare = TRUE, method = c("Joint",
   "AD"), silent = FALSE, hess = FALSE, parallel = FALSE, ...)
{
   method <- match.arg(method)
	if (pool){
			tv<- test.var.het(y)
			pv<- round(tv$p.value, 0)
			wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
			if(pv <= 0.05)	warning(wm)
	}


   if (ng == 1) {
      warning("Fitting stasis model (because ng=1)")
      if (method == "AD")
         ww <- opt.Stasis(y, pool = pool, hess = hess, ...)
      else if (method == "Joint")
         ww <- opt.joint.Stasis(y, pool = pool, hess = hess, ...)
      return(ww)
   }
   ns <- length(y$mm)
   GG <- shifts(ns, ng, minb = minb)
   nc <- ncol(GG)
   if (!silent)
      cat("Total # hypotheses: ", nc, "\n")
   i<- 1  # done to avoid R Check error in response to foreach iterators


   if (parallel==TRUE) {
      # set up cluster
      cores<-parallel::detectCores()
      cl<-parallel::makeCluster(cores-1)
      doParallel::registerDoParallel(cl)
      # run loop
      wl<- foreach::foreach (i = 1:nc, .packages=c('paleoTS')) %dopar% {
         #if (!silent)  cat(i, " ")
         gg <- shift2gg(GG[, i], ns)
         if (method == "AD")
            w <- opt.punc(y, gg, oshare = oshare, pool = pool,
               hess = hess, ...)
         else if (method == "Joint")
            w <- opt.joint.punc(y, gg, oshare = oshare, pool = pool,
               hess = hess, ...)
         w  #return w to list wl
      }
      parallel::stopCluster(cl)	      # kill cluster

   # non-parallel version
   }else if (parallel==FALSE) {
      wl<- foreach::foreach (i = 1:nc, .packages=c('paleoTS')) %do% {
         if (!silent) cat(i, " ")
         gg <- shift2gg(GG[, i], ns)
         if (method == "AD")
            w <- opt.punc(y, gg, oshare = oshare, pool = pool,
               hess = hess, ...)
         else if (method == "Joint")
            w <- opt.joint.punc(y, gg, oshare = oshare, pool = pool,
               hess = hess, ...)
         w  #return w to list wl
      }
   }

   logl<-sapply(wl, function(x) x$logL) #extract logl values from wl
   if(!silent) cat("\n")
   winner <- which.max(logl)
   wt <- wl[[winner]]

   # need to reformulate paleoTS object with altered K bc shifts now count as free parameters
   ww <- as.paleoTSfit(logL = wt$logL, parameters = wt$parameters, modelName = wt$modelName,
                       method = wt$method, K = wt$K + ng - 1, n = wt$n, se = wt$se)
   ss <- GG[, winner]
   names(ss) <- paste("shift", 1:(ng - 1), sep = "")
   ww$parameters <- append(ww$parameters, ss)
   ww$all.logl <- logl
   ww$GG <- GG
   return(ww)
}

##### GRW shift functions #####


#' Simulate (general) random walk with shift(s) in generating parameters
#'
#' @param ns vector of the number of samples in each segment
#' @param ms vector of mean step parameter in each segment
#' @param vs vector of step variance parameter in each segment
#' @param nn vector of sample sizes, one for each population
#' @param tt vector of samples times (ages)
#' @param vp phenotypic variance in each sample
#'
#' @details Simulates under a model in which a sequence is divided into two or more segments.
#' Trait evolution proceeds as a general random walk, with each segment getting its own
#' generating parameters (\code{mstep}, \code{vstep}).
#'
#' @return a \code{paleoTS} object with the simulated time-series
#' @seealso \code{\link{sim.GRW}}, \code{\link{sim.sgs}}, \code{\link{opt.GRW.shift}}
#' @export
#'
#' @examples
#' x <- sim.GRW.shift(ns = c(10,10,10), ms = c(0, 1, 0), vs = c(0.1,0.1,0.1))
#' plot(x)
#' abline(v = c(9.5, 19.5), lty = 3, lwd = 2, col = "blue")  # shows where dynamics shift
#' text (c(5, 15, 25), c(2,2,2), paste("segement", 1:3, sep =" "), col = "blue")
sim.GRW.shift <- function (ns=c(10,10), ms=c(0,1), vs=c(0.5,0.5), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
{
  nr<- length(ms)
  cns<- cumsum(ns)
  xl<- list()
  for (i in 1:nr)
   {
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else
   	 {
   	  start.i<- cns[i-1]+1
   	  end.i<- cns[i-1]+ns[i]
   	 }

   	 xl[[i]]<- sim.GRW(ns=ns[i], ms=ms[i], vs=vs[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   	 if (i>1)
   	 	xl[[i]]$mm<- xl[[i]]$mm + xl[[i-1]]$mm[length(xl[[i-1]]$mm)]
   }

  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.GRW.shift()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(ms, vs, shft)
  names(y$genpars)<- c(paste("mstep",1:nr,sep=""), paste("vstep",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
  return (y)
}

#' Fit random walk model with shift(s) in generating parameters
#'
#' @param y a \code{paloeTS} object
#' @param ng number of segments in the sequence
#' @param minb minimum number of populations in each segment
#' @param model numeric, specifies exact evolutionary model; see Details
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param silent logical, if TRUE, progress updates are suppressed
#'
#' @details Fits a model in which a sequence is divided into two or more segments and
#' trait evolution proceeds as a general random walk, with each segment (potentially)
#' getting its own generating parameters (\code{mstep}, \code{vstep}). \cr \cr
#'
#' This function tests for shifts after each population, subject to the
#' constraint that the number of populations in each segment is always >= \code{minb}. The
#' shiftpoint yielding the highest log-likelihood is returned as the solution, along with
#' the log-likelihoods (\code{all.logl} of all tested shift points (\code{GG}). \cr \cr

#'
#' Different variants of the model can be specified by the \code{model} argument:
#' \itemize{
#' \item \code{model = 1:  } \code{mstep} is separate across segments; \code{vstep} is shared
#' \item \code{model = 2:  } \code{mstep} is shared across segments; \code{vstep} is separate
#' \item \code{model = 3:  } \code{mstep} is set to zero (unbiased random walk); \code{vstep}
#' is separate across segments
#' \item \code{model = 4:  } \code{mstep} and \code{vstep} are both separate across segments
#' }
#'
#'
#' @return a \code{paleoTSfit} object
#' @seealso \code{\link{sim.GRW.shift}}
#' @export
#'
#' @examples
#' x <- sim.GRW.shift(ns = c(15,15), ms = c(0, 1), vs = c(0.1,0.1))
#' w.sep <- opt.GRW.shift(x, ng = 2, model = 4)
#' w.sameVs <- opt.GRW.shift(x, ng = 2, model = 1)
#' compareModels(w.sep, w.sameVs)
#' plot(x)
#' abline(v = x$tt[16], lwd = 3)  # actual shift point
#' abline(v = x$tt[w.sameVs$par["shift1"]], lty = 3, col = "red", lwd = 2) # inferred shift point
opt.GRW.shift<- function(y, ng=2, minb=7, model=1, pool=TRUE, silent=FALSE)
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

# internal function, not exported
opt.RW.SameMs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
  # estimates shared Ms model across multiple sequences
{
  if (class(yl)=="paleoTS")
    stop("Function opt.SameMs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)

  # pool variances if needed
  if (pool){
    for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
  {
    gg<- mle.GRW(yl[[i]])
    p0m[i]<- gg[1]
    if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
    p0v[i]<- gg[2]
  }
  p0<- c(stats::median(p0m), p0v)  # shared Ms, followed by separate Vs for each sequence
  names(p0)<- c("mstep", paste("vstep", 1:nseq, sep=""))
  if (is.null(cl$ndeps))	cl$ndeps<- rep(1e-8, length(p0))

  # optimize logL
  ll<- c(NA, rep(0,nseq))
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, yl=yl), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, yl=yl), silent=TRUE)

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameMs.Mult", method='AD', K=nseq+1, n=n, se=w$se)

  return (wc)
}

# internal function, not exported
opt.RW.SameVs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
  # estimates shared Vs model across multiple sequences
{
  if (class(yl)=="paleoTS")
    stop("Function opt.SameVs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)

  # pool variances if needed
  if (pool){
    for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
  {
    gg<- mle.GRW(yl[[i]])
    p0m[i]<- gg[1]
    if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
    p0v[i]<- gg[2]
  }
  p0<- c(p0m, stats::median(p0v))  # separate Ms, followed by shared Vs for each sequence
  names(p0)<- c(paste("mstep", 1:nseq, sep=""), "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- rep(1e-8, length(p0))

  # optimize logL
  ll<- c(rep(NA,nseq), 0)
  if (meth=="L-BFGS-B")
    w<- optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, yl=yl)
  else
    w<- optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, yl=yl)

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameVs.Mult", method='AD', K=nseq+1, n=n, se=w$se)

  return (wc)
}

# internal function, not exported
logL.SameMs <- function (p, yl)
  # computes logL over >1 sequence, of model in which all sequences have the
  # same directionality (Mstep), with different step variances
  # y is list of nseq paleoTS objects, p is array of K+1 parameters {m, v-1,..v-nseq}
{
  Sm<-0
  nseq<- length(yl)
  for (i in 1:nseq)
    Sm<- Sm + logL.GRW(p=c(p[1],p[i+1]), yl[[i]])

  return(Sm)
}

# internal function, not exported
logL.SameVs <- function (p, yl)
  # computes logL over >1 sequence, of model in which all sequences have the
  # same step variance, with different steo means
  # y is list of nseq paleoTS objects, p is array of K+1 parameters {m-1,..m-nseq, vs}
{
  Sm<-0
  nseq<- length(yl)
  for (i in 1:nseq)
    Sm<- Sm + logL.GRW(p=c(p[i],p[nseq+1]), yl[[i]])

  return(Sm)
}




#' Simulate protracted punctuation
#'
#' This function simulates a punctuated change that is is protracted enough that
#' it is captured by multiple transitional populations.  Trait evolution starts in stasis,
#' shifts to a general random walk, and then shifts back into stasis.
#'
#' @param ns vector with the number of samples in each segment
#' @param theta trait mean for initial stasis segment
#' @param omega trait variance for stasis segments
#' @param ms step mean during random walk segment
#' @param vs step variance during random walk segment
#' @param nn vector of sample sizes for each population
#' @param tt vector of times (ages) for each population
#' @param vp phenotypic trait variance for each population
#'
#' @details Trait evolution proceeds in three segments: Stasis, General random walk, stasis (sgs).
#' The initial stasis segment has a mean of \code{theta} and variance \code{omega} before
#' shifting in the second segment to a general random walk with parameters \code{ms} and
#' \code{vs}. Finally, the third segment is a return to stasis, centered around the trait value
#' of the last population of the random walk.
#'
#' @return a \code{paleoTS} object
#' @export
#'
#' @examples
#' x <- sim.sgs()  # default values OK
#' plot(x)
sim.sgs <- function (ns=c(20,20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate stasis-grw-stasis sequence, take theta2 to be final value after grw part
{
  xl<- list()
  for (i in 1:3)
   {
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1
   	 }

   	 if (i==2)
   	 	xl[[i]]<- sim.GRW(ns[2], ms, vs, nn=nn[start.i:end.i], tt=tt[start.i:end.i], vp=vp)

   	 else
   	 	xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta, omega=omega, vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }

  ## add offsets
  xl[[2]]$mm<- xl[[2]]$mm + xl[[1]]$MM[ns[1]]
  xl[[3]]$mm<- xl[[3]]$mm + xl[[2]]$MM[ns[2]]

  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.sgs()"
  y$genpars <- c(theta, omega, ms, vs)
  names(y$genpars)<- c("theta","omega", "ms","vs")
  return(y)
}


# internal function, not exported
logL.sgs<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1]
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5:6] }
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4:5]	}

  l1<- logL.Stasis(p=c(th[1], om[1]), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om[2]), yl[[3]])

  logl<- l1+l2+l3
  return(logl)
}

# internal function, not exported
logL.sgs.omega<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model, omega shared over stasis segments
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1]
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5]	}
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4]	}

  l1<- logL.Stasis(p=c(th[1], om), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om), yl[[3]])

  logl<- l1+l2+l3
  return(logl)
}

# internal function, not exported
opt.sgs<- function(y,gg,cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare=TRUE, model="GRW")
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
 names(p0)<- pn

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
  else 			w$se<- NULL

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('SGS', model, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)
  return(wc)
}



#' Fit a model of trait evolution with a protracted punctuation.
#'
#' This function fits a model of punctuated change that is is protracted enough that
#' it is captured by multiple transitional populations.  Trait evolution starts in stasis,
#' shifts to a general random walk, and then shifts back into stasis.
#'
#' @param y a \code{paleoTS} object
#' @param minb minimum number of populations within each segment
#' @param oshare logical, if TRUE, variance assumed to be shared (equal) across segments
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param silent logical, if TRUE, progress updates are suppressed
#' @param hess if TRUE, standard errors computed from the Hessian matrix are returned
#' @param meth optimization method, passes to \code{optim}
#' @param model type of random walk: \code{"URW"}, unbiased random walk, or \code{"GRW"},
#' a general (directional) random walk
#'
#' @return  a \code{paleoTSfit} object
#' @seealso \code{\link{fitGpunc}}
#' @export
#'
#' @examples
#' x <- sim.sgs(ns = c(15, 15, 15))  # default values OK
#' w <- fit.sgs(x, minb = 10)  # increase minb so example takes less time; not recommended!
#' plot(x)
#' abline(v = c(16, 31), lwd = 3)  # actual shifts
#' abline(v = c(w$parameters[6:7]), lwd = 2, lty = 3, col = "red")  # inferred shifts

fit.sgs<- function(y, minb=7, oshare=TRUE, pool=TRUE, silent=FALSE, hess=FALSE, meth="L-BFGS-B", model="GRW")
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




##### Complex mode functions from PNAS #####

# internal function, not exported
logL.joint.URW.Stasis<- function(p, y, gg)
 {
 	# preliminaries
 	anc<- p[1]
 	vs<- p[2]
 	theta<- p[3]
 	omega<- p[4]
 	n<- length(y$mm)

	# get vector of means
	st.seg<- which(gg==2)
 	rw.seg<- which(gg==1)
 	tt.rw<- y$tt[rw.seg]
	M<- c(rep(anc, length(rw.seg)), rep(theta, length(st.seg)))
	M<- unname(M)

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<-  VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")

	# logL from mvn
	S<- mnormt::dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }


# internal function, not exported
 logL.joint.GRW.Stasis<- function(p, y, gg)
 {
 	# preliminaries
 	anc<- p[1]
 	ms<- p[2]
 	vs<- p[3]
 	theta<- p[4]
 	omega<- p[5]
 	n<- length(y$mm)

	# get vector of means
 	#ts<- max(y$tt[gg==1])  # time of shift point
	st.seg<- which(gg==2)
 	rw.seg<- which(gg==1)
 	tt.rw<- y$tt[rw.seg]
	M<- c(rep(anc, sum(gg==1)) + ms*y$tt[gg==1], rep(theta, length(st.seg)))
	M<- unname(M)

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<-  VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")


	# logL from mvn
	S<- mnormt::dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }



 # internal function, not exported
 opt.joint.RW.Stasis<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)

  # get initial parameter estimates
  small<- 1e-8
  if(rw.model=="URW")  {p0rw<- mle.URW(sub.paleoTS(y, ok=gg==1)); K<- 5} # assumes shift point is free parameter
  else				   {p0rw<- mle.GRW(sub.paleoTS(y, ok=gg==1)); K<- 6}
  p0st<- mle.Stasis(sub.paleoTS(y, ok=gg==2))
  if(p0rw["vstep"] <= small)	p0rw["vstep"]<- 100*small
  if(p0st["omega"] <= small) 	p0st["omega"]<- 100*small
  p0anc<- y$mm[1]
  names(p0anc)<- "anc"
  p0<- c(p0anc, p0rw, p0st)
  #print(p0)

  #cl$parscale <-
  ll.urw<- c(NA,small,NA,small)
  ll.grw<- c(NA,NA,small,NA,small)
  if(rw.model=="URW")	{ ll<- ll.urw}
  else					{ ll<- ll.grw}
  if(meth!="L-BFGS-B")	ll<- NULL  # sets so to determine meth

  if(rw.model=="URW")	{
  		w <- try(optim(p0, fn=logL.joint.URW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,100,1,10))
  				w <- try(optim(p0, fn=logL.joint.URW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
  				}
  } else if(rw.model=="GRW"){
  		w <- try(optim(p0, fn=logL.joint.GRW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y) , silent=TRUE)
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,100,1,10))
 				w <- try(optim(p0, fn=logL.joint.GRW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y) , silent=TRUE)
 				}}
  # add more information to results
  if(class(w)=="try-error")	{
  		wc<- as.paleoTSfit(logL=NA, parameters=NA, modelName=paste(rw.model, "Stasis", sep='-'), method='Joint', K=K, n=length(y$mm), se=NULL)
  		return(wc)
  		}
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(rw.model, "Stasis", sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)
 }

 # internal function, not exported
 opt.AD.RW.Stasis<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)

  # split and optimize AD separately for segments
  yl<- split4punc(y, gg)
  if(rw.model=="URW")	{ w1<- opt.URW(yl[[1]], pool=pool, meth=meth, hess=hess); K<- 4}
  else 					{ w1<- opt.GRW(yl[[1]], pool=pool, meth=meth, hess=hess); K<- 5}
  w2<- opt.Stasis(yl[[2]], pool=pool, meth=meth, hess=hess)


   # add more information to results
  if (hess)	{	se1<- sqrt(diag(-1*solve(w1$hessian)))
  				se2<- sqrt(diag(-1*solve(w2$hessian)))
  				se<- c(se1, se2)  }
  else			se<- NULL
  wc<- as.paleoTSfit(logL=w1$logL+w2$logL, parameters=c(w1$par, w2$par), modelName=paste(rw.model, "Stasis", sep='-'), method='AD', K=K, n=length(y$mm)-1, se=se)

  return(wc)
 }

 # internal function, not exported
  logL.joint.Stasis.URW<- function(p, y, gg)
 {
 	# preliminaries
 	#cat(round(p,6), "\t[")
 	theta<- p[1]
 	omega<- p[2]
 	vs<- p[3]
 	n<- length(y$mm)

	# get vector of means
	M<- rep(theta, n)
	M<- unname(M)

 	# compute covariance matrix
 	st.seg<- which(gg==1)
 	VVst<- diag(omega, nrow=length(st.seg))
 	rw.seg<- which(gg==2)
 	tt.rw<- y$tt[rw.seg] - y$tt[rw.seg[1]-1]
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<- VVst
 	VVtot[rw.seg, rw.seg]<- VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error

	# logL from mvn
	S<- mnormt::dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }


  # internal function, not exported
 logL.joint.Stasis.GRW<- function(p, y, gg)
 {
 	# preliminaries
 	theta<- p[1]
 	omega<- p[2]
 	ms<- p[3]
 	vs<- p[4]
 	n<- length(y$mm)

	# get vector of means
 	st.seg<- which(gg==1)
 	rw.seg<- which(gg==2)
 	#ts<- max(y$tt[gg==1])  # time of shift point
 	tt.rw<- y$tt[rw.seg] - y$tt[rw.seg[1]-1]
	M<- c(rep(theta, length(st.seg)), theta + ms*tt.rw)
	M<- unname(M)

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<- VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")

	# logL from mvn
	S<- mnormt::dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }


 # internal function, not exported
 opt.joint.Stasis.RW<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)

  # get initial parameter estimates
  small<- 1e-8
  if(rw.model=="URW")  {p0rw<- mle.URW(sub.paleoTS(y, ok=gg==2)); K<- 4} # assumes shift point is free parameter
  else				   {p0rw<- mle.GRW(sub.paleoTS(y, ok=gg==2)); K<- 5}
  p0st<- mle.Stasis(sub.paleoTS(y, ok=gg==1))
  if(p0rw["vstep"] <= small)	p0rw["vstep"]<- 100*small
  if(p0st["omega"] <= small) 	p0st["omega"]<- 100*small
  p0<- c(p0st, p0rw)
  #cat(p0, "\n\n")


  ll.urw<- c(NA,small,small)
  ll.grw<- c(NA,small,NA,small)
  if(rw.model=="URW")	{ ll<- ll.urw}
  else					{ ll<- ll.grw}
  if(meth!="L-BFGS-B")	ll<- -Inf  # sets so to determine meth


  if(rw.model=="URW")	{
  		w <- try(optim(p0, fn=logL.joint.Stasis.URW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,100))
		  		w <- try(optim(p0, fn=logL.joint.Stasis.URW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
  				}
  }else if(rw.model=="GRW")	{
  		w <- try(optim(p0, fn=logL.joint.Stasis.GRW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,10,100))
		  		w <- try(optim(p0, fn=logL.joint.Stasis.GRW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)
 				}}

  # add more information to results
  if(class(w)=="try-error")	{
  		wc<- as.paleoTSfit(logL=NA, parameters=NA, modelName=paste("Stasis", rw.model, sep='-'), method="Joint", K=K, n=length(y$mm), se=NULL)
  		return(wc)
  		}

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste("Stasis", rw.model, sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)
 }

 # internal function, not exported
 opt.AD.Stasis.RW<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)

  # split and optimize AD separately for segments
  yl<- split4punc(y, gg)
  w1<- opt.Stasis(yl[[1]], pool=pool, meth=meth, hess=hess)
  if(rw.model=="URW")	{ w2<- opt.URW(yl[[2]], pool=pool, meth=meth, hess=hess); K<- 4}
  else 					{ w2<- opt.GRW(yl[[2]], pool=pool, meth=meth, hess=hess); K<- 5}



   # add more information to results
  if (hess)	{	se1<- sqrt(diag(-1*solve(w1$hessian)))
  				se2<- sqrt(diag(-1*solve(w2$hessian)))
  				se<- c(se1, se2)  }
  else			se<- NULL
  wc<- as.paleoTSfit(logL=w1$logL+w2$logL, parameters=c(w1$par, w2$par), modelName=paste("Stasis", rw.model, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=se)

  return(wc)
 }


#' Simulate trait evolution with a mode shift
#'
#' Trait evolution is modeled as a shift from a random walk (general or unbiased)
#' to stasis, or vice versa.
#'
#' @param ns vector of the number of samples in each segment
#' @param order whether stasis or random walk come first, one of \code{"Stasis-RW"} or
#'   \code{"RW-Stasis"}
#' @param anc starting trait value
#' @param omega variance of stasis segment
#' @param ms step mean during random walk segment
#' @param vs step variance during random walk segment
#' @param vp phenotypic trait variance for each population
#' @param nn vector of sample sizes for each population
#' @param tt vector of times (ages) for each population
#'
#' @details The \code{anc} argument is the starting trait value, and if the
#' first segment is stasis, this is also the value of the stasis mean. When the first segment
#' is a random walk, the stasis mean in the second segment is equal to the true trait mean at
#' the end of the initial random walk.
#'
#' @return a \code{paleoTSfit} object
#' @seealso \code{\link{fitModeShift}}
#' @export
#'
#' @examples
#' x1 <- sim.Stasis.RW(omega = 0.1, ms = 5, order = "Stasis-RW")
#' x2 <- sim.Stasis.RW(omega = 0.1, ms = 5, order = "RW-Stasis")
#' plot(x1)
#' plot(x2, add = TRUE, col = "blue")
#' abline(v = 19, lty=3)
#'
 sim.Stasis.RW<- function(ns=c(20,20), order = c("Stasis-RW", "RW-Stasis"), anc=0, omega=1, ms=0, vs=1, vp=1, nn=30, tt=NULL)
 {
   # Add a sample to the second segments because its first sample will be deleted
   # Set-up tt for sequences accordingly
   ns[2] <- ns[2] + 1
   ttl<- list()
   if(is.null(tt)) {
     ttl[[1]]<- 0:(ns[1]-1)
     ttl[[2]]<- (ns[1]-1):(sum(ns)-2)
   } else {
     ttl[[1]] <- tt[1:(ns[1])]
     ttl[[2]] <- tt[(ns[1]+1):sum(ns)]
     }


   # check with model comes first
   order = match.arg(order)
   if (order == "Stasis-RW"){
     st <- sim.Stasis(ns = ns[1], theta = anc, omega = omega, vp=vp, nn=rep(nn, ns[1]), tt=ttl[[1]])
     rw <- sim.GRW(ns=ns[2], ms=ms, vs=vs, vp=vp, nn=rep(nn, ns[2]), tt=ttl[[2]])

     # shift rw so trait values start at ending value of st; remove first population of rw
     rw$MM <- rw$MM + st$MM[ns[1]]
     rw$mm <- rw$mm + st$MM[ns[1]]
     rw <- sub.paleoTS(rw, ok = 2:ns[2], reset.time = FALSE)

     y <- cat.paleoTS(list(st, rw))
     } else {
     rw <- sim.GRW(ns=ns[1], ms=ms, vs=vs, vp=vp, nn=rep(nn, ns[1]), tt=ttl[[1]])
     st <- sim.Stasis(ns = ns[2], theta = rw$MM[ns[1]], omega = omega, vp=vp, nn=rep(nn, ns[2]), tt=ttl[[2]])

     # st already starts at end point of rw; remove first population of st
     st <- sub.paleoTS(st, ok = 2:ns[2], reset.time = FALSE)

     y <- cat.paleoTS(list(rw, st))
    }

   # now drop 1st point of s2 and concatenate
   y$genpars<- c(anc, omega, ms, vs, ns[1]+1)
   names(y$genpars)<- c("anc", "omega", "ms", "vs", "shift.start")
   y$label<- "Created by sim.Stasis.RW()"

   return(y)
 }


#' Fit model in which the mode of trait evolution shifts once
#'
#' Trait evolution is modeled as a shift from a random walk (general or unbiased) to stasis, or
#' vice versa.
#'
#' @param y  \code{paleoTS} object
#' @param minb minimum number of populations within each segment
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param order whether stasis or random walk come first, one of \code{Stasis-RW} or
#'   \code{RW-Stasis}
#' @param rw.model whether the random walk segment is an unbiased random walk, \code{URW}
#'   or a general random walk, \code{GRW}
#' @param method parameterization to use: \code{Joint} or \code{AD}
#' @param silent logical, if TRUE, progress updates are suppressed
#' @param hess if TRUE, standard errors computed from the Hessian matrix are returned
#' @param ... other arguments, passed to optimization functions
#'
#' @return a \code{paleoTSfit} object
#' @seealso \code{\link{sim.Stasis.RW}}
#' @export
#'
#' @examples
#' x <- sim.Stasis.RW(ns = c(15, 15), omega = 0.5, ms = 1, order = "Stasis-RW")
#' plot(x)
#' w <- fitModeShift(x, order = "Stasis-RW", rw.model = "GRW")
#' abline(v = x$tt[15], lwd = 3)  # actual shift point
#' abline(v = x$tt[w$par["shift1"]], lty = 3, lwd = 2, col = "red") # inferred shift
#'
 fitModeShift<- function(y, minb=7, pool=TRUE, order=c("Stasis-RW", "RW-Stasis"), rw.model=c("URW","GRW"), method=c('Joint', 'AD'), silent=FALSE, hess=FALSE,...)
## optimize models (with some min n per segment) that shift from GRW/URW to Stasis, or vice versa
{
 method<- match.arg(method)
 order<- match.arg(order)
 rw.model<- match.arg(rw.model)

 ns<- length(y$mm)
 ng<- 2
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("Total # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)

 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    # the different gg for AD and J is required for the interpretation of the "shift" parameters to be the same across parameterizations
    gg<- shift2gg(GG[,i], ns)

    if(method=='AD') {
    		if(order=="Stasis-RW")  		w<- opt.AD.Stasis.RW(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    		else if (order=="RW-Stasis")	w<- opt.AD.RW.Stasis(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    	}
    if (method=='Joint'){
    		if(order=="Stasis-RW")  		w<- opt.joint.Stasis.RW(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    		else if (order=="RW-Stasis")	w<- opt.joint.RW.Stasis(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
		}

    logl[i]<- w$logL
    wl[[i]]<- w
  }
 if(!silent)	cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}


#' Bootstrap test to see if a complex model is significantly better than a simple
#' model.
#'
#'
#' @param y a \code{paleoTS} object
#' @param simpleFit a \code{paleoTSfit} object, representing the model fit of a
#'   simple model
#' @param complexFit a \code{paleoTSfit} object, representing the model fit of a
#'   complex model
#' @param nboot number of replications for parametric bootstrapping
#' @param minb minimum number of populations within each segment
#' @param ret.full.distribution logical, indicating if the null distribution for
#'   the likelihood ratio from the parametric bootstrap should be returned
#' @param parallel logical, if TRUE, the bootstrapping is done using parallel
#'   computing
#' @param ... further arguments, passed to optimization functions
#'
#' @details Simulations suggest that AICc can be overly liberal with complex
#' models with mode shifts or punctuations (Hunt et al., 2015). This function
#' implements an alternative of parametric boostrapping to compare the fit of a
#' simple model with a complex model. It proceeds in five steps: \enumerate{
#' \item Compute the observed gain in support from the simple to complex model
#' as the likelihood ratio, \eqn{LR_obs = -2(logL_simple - logL_complex) } \item
#' Simulate trait evolution under the specified simple model \code{nboot} times
#' \item Fit to each simulated sequence the specified simple and complex models
#' \item Measure the gain in support from simple to complex as the bootstrap
#' likelihood ratio for each simulated sequence \item Compute the P-value as the
#' percentile of the bootstrap distribution corresponding to the observed LR. }
#'
#' Argument \code{simpleFit} should be a \code{paleoTS} object returned by the
#' function \code{fitSimple} or similar functions (e.g., \code{opt.joint.GRW,
#' opt.GRW}, etc.). Argument \code{complexFit} must be a \code{paleoTS} object
#' returned by \code{fitGpunc} or \code{fitModeShift}.
#'
#' Calculations can be speeded up by setting \code{parallel = TRUE}, which uses
#' functions from the \code{\link{doParallel}} package to run the bootstrap
#' replicates in parallel, using one fewer than the number of detected cores.
#'
#' @return A list of the observed likelihood ratio statistic, \code{LRobs}, the
#'   P-value of the test, and the number of bootstrap replicates. If
#'   \code{ret.full.distribution = TRUE}, the null distribution of likelihood
#'   ratios generated by parametric bootstrapping is also returned.
#' @seealso \code{\link{sim.Stasis.RW}}, \code{\link{fitModeShift}}
#' @references Hunt, G., M. J. Hopkins and S. Lidgard. 2015. Simple versus
#' complex models of trait evolution and stasis as a response to environmental
#' change. PNAS 112(16): 4885-4890.

#' @export
#'
#' @examples
#' \dontrun{
#' x <- sim.Stasis.RW(ns = c(15, 15), omega = 0.5, ms = 1, order = "Stasis-RW")
#' ws <- fitSimple(x)
#' wc <- fitModeShift(x, order = "Stasis-RW", rw.model = "GRW")
#' bootSimpleComplex(x, ws, wc, nboot = 50, minb = 7)  # nboot too low for real analysis!
#' }
bootSimpleComplex<- function(y, simpleFit, complexFit, nboot=99, minb=7, ret.full.distribution=FALSE, parallel=FALSE, ...)
# function to to parametric bootstrapping to test if complex model is better than a simple model
# AICc seems often too liberal in favoring complex models
# consider comparing simpleFit to fit of all complex models?
{
	# make sure models to be tested are available
	sName<- simpleFit$modelName
	cName<- complexFit$modelName
	sAvailModels<- c("StrictStasis", "Stasis", "URW", "GRW")
	cAvailModels<- c("Punc-1", "URW-Stasis", "GRW-Stasis", "Stasis-URW", "Stasis-GRW")
	if (! (sName %in% sAvailModels) )	stop(paste("Simple model [", sName, "]", "not among available models [", sAvailModels, "]"))
	if (! (cName %in% cAvailModels) )	stop(paste("Complex model [", cName, "]", "not among available models [", cAvailModels, "]"))


	# compute observed LR
	LR<- -2*(simpleFit$logL - complexFit$logL)

	# generate bootstrap samples, fit simple and complex
	xboot<- list()
	pp<- simpleFit$par
	ns<- length(y$mm)
	vp<- pool.var(y)

   if (parallel==TRUE) {
      cores<-parallel::detectCores()
      cl<-parallel::makeCluster(cores-1)
      doParallel::registerDoParallel(cl)

      wl<- foreach::foreach (i = 1:nboot, .packages=c('paleoTS')) %dopar% {

		# create simulated dataset, fit simple model
		if (sName == "GRW"){
				xboot<- sim.GRW(ns=ns, ms=pp['mstep'], vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.GRW(xboot, ...)	}
		if (sName == "URW") {
				xboot<- sim.GRW(ns=ns, ms=0, vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.URW(xboot, ...)	}
		if (sName == "Stasis"){
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=pp['omega'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.Stasis(xboot, ...)	}
		if (sName == "StrictStasis") {
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=0, vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.StrictStasis(xboot, ...)	}

		# fit complex model
		if(cName=="Punc-1")	compBFit<- fitGpunc(xboot, ng=2, method="Joint", silent=TRUE, ...)
		if(cName=="URW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="GRW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="GRW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-URW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-GRW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="GRW", method="Joint", silent=TRUE, ...)

		bres<- list(simpBFit=simpBFit, compBFit=compBFit, xboot=xboot)
		bres
      }
      parallel::stopCluster(cl)	      # kill cluster
	}
	if(parallel==FALSE){
	      wl<- foreach::foreach (i = 1:nboot, .packages=c('paleoTS')) %do% {

		# create simulated dataset, fit simple model
		if (sName == "GRW"){
				xboot<- sim.GRW(ns=ns, ms=pp['mstep'], vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.GRW(xboot, ...)	}
		if (sName == "URW") {
				xboot<- sim.GRW(ns=ns, ms=0, vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.URW(xboot, ...)	}
		if (sName == "Stasis"){
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=pp['omega'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.Stasis(xboot, ...)	}
		if (sName == "StrictStasis") {
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=0, vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.StrictStasis(xboot, ...)	}

		# fit complex model
		if(cName=="Punc-1")	compBFit<- fitGpunc(xboot, ng=2, method="Joint", silent=TRUE, ...)
		if(cName=="URW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="GRW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="GRW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-URW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-GRW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="GRW", method="Joint", silent=TRUE, ...)

		bres<- list(simpBFit=simpBFit, compBFit=compBFit, xboot=xboot)
		bres
      }
	}

	#return(wl)

	# compile information
	logL.simp<- sapply(wl, FUN=function(x) x$simpBFit$logL)
	logL.comp<- sapply(wl, FUN=function(x) x$compBFit$logL)
	LRnull<- -2* (logL.simp - logL.comp)
	cc<- is.finite(LRnull)
	p.value<- (sum(LRnull[cc] > LR) +1)  / (length(LRnull[cc])+1)
	res<- list(LRobs=LR, p.value=p.value, nboot = nboot)

	if(ret.full.distribution){
			res$nullLR<- LRnull[cc]	}
			#res$boot<- 	lapply(wl, FUN=function(x) x$xboot)	}

	return(res)
}


#' Fit large set of models to a time-series
#'
#' This function fits nine models to a time-series following Hunt et al. (2015). These
#' include the simple models fit by \code{fit4models} along with mode shift and
#' punctuation models.
#'
#' @param y a \code{paleoTS} object
#' @param silent logical, if TRUE, progress updates are suppressed
#' @param method parameterization to use: \code{Joint} or \code{AD}; see Details
#' @param ... other arguments, passed to optimization functions
#'
#' @references
#' Hunt, G., M. J. Hopkins and S. Lidgard. 2015. Simple versus complex models of trait evolution
#' and stasis as a response to environmental change. PNAS 112(16): 4885-4890.

#' @return if silent = FALSE, a table of model fit statistics, also printed to the
#' screen.  if silent = TRUE, a list of the model fit statistics and model parameter values.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- sim.Stasis.RW(ns = c(15, 15), omega = 0.5, ms = 1, order = "Stasis-RW")
#' plot(x)
#' fit9models(x)
#' }
fit9models<- function(y, silent=FALSE, method=c("Joint", "AD"), ...)
{
  args<- list(...)
  check.var<- TRUE
  if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
  if (check.var){
    tv<- test.var.het(y)
    pv<- round(tv$p.value, 0)
    wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
    if(pv <= 0.05)	warning(wm)
  }

  method <- match.arg(method)
  if (method == "AD") {
    if(!silent) cat("Fitting simple models...\n")
    m4 <- opt.GRW(y, ...)
    m3 <- opt.URW(y, ...)
    m2 <- opt.Stasis(y, ...)
    m1 <- opt.StrictStasis(y, ...)
  }
  else if (method == "Joint") {
    if(!silent)	cat("Fitting simple models...\n")
    m4 <- opt.joint.GRW(y, ...)
    m3 <- opt.joint.URW(y, ...)
    m2 <- opt.joint.Stasis(y, ...)
    m1 <- opt.joint.StrictStasis(y, ...)
  }

  if(!silent)	cat("Fitting punctuational model...\n")
  m5 <- fitGpunc(y, ng=2, method=method, silent=silent, ...)
  if(!silent)	cat("Fitting Stasis-URW model...\n")
  m6 <- fitModeShift(y, order="Stasis-RW", rw.model="URW", method=method, silent=silent, ...)
  if(!silent)	cat("Fitting Stasis-GRW model...\n")
  m7 <- fitModeShift(y, order="Stasis-RW", rw.model="GRW", method=method, silent=silent, ...)
  if(!silent)  cat("Fitting URW-Stasis model...\n")
  m8 <- fitModeShift(y, order="RW-Stasis", rw.model="URW", method=method, silent=silent, ...)
  if(!silent)  cat("Fitting GRW-Stasis model...\n")
  m9 <- fitModeShift(y, order="RW-Stasis", rw.model="GRW", method=method, silent=silent,...)


  mc <- compareModels(m1, m2, m3, m4, m5, m6, m7, m8, m9, silent = silent)
  invisible(mc)
}






