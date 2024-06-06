# multModels #

#' Fit the same simple model across multiple time-series
#'
#' @param yl a list of \code{paleoTS} objects
#' @param model the model to fit; see Details
#' @param method parameterization to use: \code{Joint} or \code{AD}
#' @param pool if TRUE, sample variances are substituted with their pooled
#'   estimate
#' @param zl for the \code{covTrack} model only, a list of covariate vectors,
#'   one each \code{paleoTS} object in \code{yl}
#' @param hess if TRUE, standard errors computed from the Hessian matrix are
#'   returned
#'
#' @details This function fits a model with shared parameters across multiple
#'   trait time-series. The most likely application would be to model a common
#'   evolutionary dynamic across different sequences, perhaps representing
#'   time-series of the same trait and lineage from different localities or time
#'   intervals.
#'
#'   Four simple models are currently implemented: \itemize{
#'   \item \strong{GRW}:
#'   parameters \code{mstep} and \code{vstep} of the general random walk are
#'   shared across sequences.
#'   \item \strong{URW}: parameter \code{vstep} of the
#'   unbiased random walk is shared across sequences.
#'   \item \strong{Stasis}:
#'   parameter \code{omega} of stasis is shared across sequences.
#'   \item
#'   \strong{covTrack}: parameters \code{b0}, \code{b1},  and \code{evar} of the
#'   covariate-tracking model are shared across sequences. }
#'
#' Under the joint parameterization, \code{method = "Joint"}, an additional parameter, \code{anc} is
#' fit, representing the ancestral (starting) trait value. This parameter is estimated separately
#' in each sequence so it is not assumed that they all start at the same trait value.
#'
#' @return a \code{paleoTSfit} object with the results of the model-fitting
#' @note The models are described in the help for \code{fitSimple} and the functions
#' linked from there.
#' @seealso \code{\link{fitSimple}}
#' @export
#'
#' @examples
#' x1 <- sim.GRW(ms = 1, vs = 0.2)
#' x2 <- sim.GRW(ms = 1, vs = 0.2)
#' fitMult(list(x1, x2), model = "GRW")
fitMult<- function(yl, model=c("GRW", "URW", "Stasis", "covTrack"), method=c("Joint", "AD"), pool=TRUE, zl=NULL, hess=FALSE)
{
  model<- match.arg(model)
  method<- match.arg(method)

  if(model=="covTrack"){
    if(method=="Joint")	w<- opt.joint.covTrack.Mult(yl, zl, pool=pool, hess=hess)
    if(method=="AD")	w<- opt.covTrack.Mult(yl, zl, pool=pool, hess=hess)
  }	else{
    if(method=="AD")		w<- opt.Mult(yl, model=model, pool=pool, hess=hess)
    if(method=="Joint")	w<- opt.joint.Mult(yl, model=model, pool=pool, hess=hess)
  }

  return(w)
}

# internal function, not exported
logL.Mult.covTrack<- function (p, yl, zl)
  # y is a _list_ of paleoTS objects, and z is a _list_ of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
  {	Smult<- Smult + logL.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)
}

# internal function, not exported
opt.covTrack.Mult<- function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
  # y and z are lists of paleoTS, and covariates, respectively
{
  if (inherits(yl,"paleoTS"))
  { stop("opt.covTrack.Mult is only for multiple paleoTS sequences\n") }
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


  regMat<- array(dim=c(nseq, 3))  # hold results of regressions for sequences separately
  colnames(regMat)<- c("b0", "b1", "evar")
  for(i in 1:nseq)
  {
    w<- stats::lm(diff(yl[[i]]$mm) ~ zl[[i]])
    aa<- stats::anova(w)
    regMat[i,]<- c(stats::coef(w), aa$Sum[2]/aa$Df[2])  # estimates of b0, b1, and evar
  }
  p0<- apply(regMat, 2, stats::median)  # initial estimates are median of separate regressions
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(yl[[1]], p0)

  w<- optim(p0, fn=logL.Mult.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,1e-6), hessian=hess, yl=yl, zl=zl)

  if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
  else		w$se<- NULL

  ff<- function(x) length(x$mm)-1
  n<- sum(sapply(yl, ff))

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='covTrack.Mult',
                     method='AD', K=2, n=n, se=w$se, convergence=w$convergence, logLFunction="logL.Mult.covTrack")
  return(wc)
}

# internal function, not exported
logL.Mult.joint.covTrack<- function (p, yl, zl)
  # y is a _list_ of paleoTS objects, and z is a *list* of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
  {	Smult<- Smult + logL.joint.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)
}

# internal function, not exported
opt.joint.covTrack.Mult<- function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
  # y and z are lists of paleoTS, and covariates, respectively
{
  if (inherits(yl, "paleoTS"))
  { stop("opt.covTrack.Mult is only for multiple paleoTS sequences\n") }
  nseq <- length(yl)
  if (pool) {
    for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE) }

  # check lengths of z
  for (i in 1:nseq)
  {
    ns<- length(yl[[i]]$mm)
    if (length(zl[[i]]) != ns)  stop("Covariate length [", length(zl[[i]]), "] does not match the sequence length [", ns, "]\n" )

  }

  regMat<- array(dim=c(nseq, 3))  # hold results of regressions for sequences separately
  colnames(regMat)<- c("b0", "b1", "evar")
  for(i in 1:nseq)
  {
    w<- stats::lm(yl[[i]]$mm ~ zl[[i]])
    aa<- stats::anova(w)
    regMat[i,]<- c(stats::coef(w), aa$Sum[2]/aa$Df[2])  # estimates of b0, b1, and evar
  }
  p0<- apply(regMat, 2, stats::median)  # initial estimates are median of separate regressions
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(yl[[1]], p0)

  w<- optim(p0, fn=logL.Mult.joint.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,NA,1e-6), hessian=hess, yl=yl, zl=zl)

  if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
  else		w$se<- NULL

  ff<- function(x) length(x$mm)
  n<- sum(sapply(yl, ff))

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='covTrack.Mult', method='Joint',
                     K=3, n=n, se=w$se, convergence=w$convergence, logLFunction="logL.Mult.joint.covTrack")
  return(wc)
}


# internal function, not exported
logL.Mult<- function (p, yl, model=c("GRW", "URW", "Stasis"))
  # calculate logL over multiple sequences
  #  here, yl is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(yl)
  for (i in 1:nseq)
  {
    if (model=="URW")
      Smult<- Smult + logL.URW(p,yl[[i]])
    else if (model=="GRW")
      Smult<- Smult + logL.GRW (p,yl[[i]])
    else if (model=="Stasis")
      Smult<- Smult + logL.Stasis(c(p[i], p[nseq]), yl[[i]])
  }
  return (Smult)
}


# internal function, not exported
opt.Mult<- function (yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
  # estimates single model across multiple sequences
  # pool=TRUE will pool variances _within_ sequences
{
  if (inherits(yl, "paleoTS"))
    stop("opt.RW.mult is only for multiple paleoTS sequences\n")
  nseq<- length(yl)

  # pool variances if needed
  if (pool){
    for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
    ll<- c(NA,1e-6)
    pp<- sapply(yl, mle.GRW)
    p0<- apply(pp, 1, stats::median)	}
  else if (model=="URW")
  { ll<- 1e-6
  pp<- sapply(yl, mle.URW)
  p0<- stats::median(pp)
  names(p0)<- "vstep"}
  else if(model=="Stasis")
  {  ll<- c(rep(NA, nseq), 1e-6)
  pp<- sapply(yl, mle.Stasis)
  p0<- c(pp[1,], stats::median(pp[2,]))
  names(p0)<- c(paste("theta", 1:nseq, sep=""), "omega")
  }
  K<- length(p0)
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(yl[[1]], p0)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE)

  # add more information to results (p0, SE, K, n, IC scores)
  nn<- sapply(yl, FUN=function(x) length(x$mm))
  n<- sum(nn) - nseq
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(model,'.Mult', sep=''),
                     method="AD", K=K, n=n, se=w$se, convergence=w$convergence, logLFunction="logL.Mult")

  return (wc)
}

# internal function, not exported
opt.joint.Mult<- function (yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
  # estimates single model across multiple sequences
  # pool=TRUE will pool variances _within_ sequences
{
  if (inherits(yl, "paleoTS") || !inherits(yl[[1]],"paleoTS"))
    stop("opt.joint.Mult() is only for a list of multiple paleoTS sequences\n")
  nseq<- length(yl)

  # pool variances if needed
  if (pool){
    for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
    ll<- c(rep(NA, nseq), NA,1e-6)
    anc0<- sapply(yl, FUN=function(x) x$mm[1])
    P0<- sapply(yl, FUN=mle.GRW)
    p0<- c(anc0, apply(P0, 1, stats::median))
    if(p0[nseq+2]<=0)	p0[nseq+2]<- 1e-5
    names(p0)<- c(paste("anc", 1:nseq, sep=""), "mstep", "vstep")
    K<- nseq+2	 	}
  else if (model=="URW"){
    ll<- c(rep(NA, nseq),1e-6)
    anc0<- sapply(yl, FUN=function(x) x$mm[1])
    P0<- sapply(yl, FUN=mle.URW)
    p0<- c(anc0, stats::median(P0))
    if(p0[nseq+1]<=0)	p0[nseq+1]<- 1e-5
    names(p0)<- c(paste("anc", 1:nseq, sep=""), "vstep")
    K<- nseq+1	 	}
  else if (model=="Stasis"){
    ll<- c(rep(NA, nseq),1e-6)
    P0<- sapply(yl, FUN=mle.Stasis)
    p0<- c(P0[1,], stats::median(P0[2,]))
    if(p0[nseq+1]<=0)	p0[nseq+1]<- 1e-5
    names(p0)<- c(paste("theta", 1:nseq, sep=""), "omega")
    K<- nseq+1	 	}

  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(yl[[1]], p0)

  if (meth=="L-BFGS-B")
    w<- optim(p0, fn=logL.joint.Mult, method=meth, lower=ll, control=cl, hessian=hess, yl=yl, model=model)
  else     w<- optim(p0, fn=logL.joint.Mult, method=meth, control=cl, hessian=hess, yl=yl, model=model)

  # add more information to results (p0, SE, K, n, IC scores)
  n<- sum(sapply(yl, FUN=function(x)length(x$mm)))

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(model,'.Mult', sep=''), method='Joint',
                     K=K, n=n, se=w$se, convergence=w$convergence, logLFunction="logL.joint.Mult")

  return (wc)
}

# internal function, not exported
logL.joint.Mult<- function (p, yl, model=c("GRW", "URW", "Stasis"))
  # calculate logL over multiple sequences
  #  here, yl is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(yl)
  for (i in 1:nseq)
  {
    if (model=="URW")
      Smult<- Smult + logL.joint.URW(p=c(p[i], p[nseq+1]), yl[[i]])	# first nseq are anc_i, then others
    else if (model=="GRW")
      Smult<- Smult + logL.joint.GRW (p=c(p[i], p[nseq+1], p[nseq+2]), yl[[i]])
    else if (model=="Stasis")
      Smult<- Smult + logL.joint.Stasis (p=c(p[i], p[nseq+1]), yl[[i]])
  }
  return (Smult)
}





