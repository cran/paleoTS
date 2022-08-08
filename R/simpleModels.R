# simpleModels #

#' Simulate random walk or directional time-series for trait evolution
#'
#' @param ns number of populations in the sequence
#' @param ms mean of evolutionary "steps"
#' @param vs variance of evolutionary "steps"
#' @param vp phenotypic variance within populations
#' @param nn vector of population sample sizes
#' @param tt vector of population times (ages)
#'
#' @return a \code{paleoTS} object
#' @details The general random walk model considers time in discrete steps.
#' At each time step, an evolutionary change is drawn at random from a distribution of
#' possible evolutionary "steps."  It turns out that the long-term dynamics of an evolving
#' lineage depend only on the mean and variance of this step distribution.  The former,
#' \code{mstep}, determined the directionality in a sequence and the latter, \code{vstep},
#' determines its volatility.
#' @note This function simulates an unbiased random walk if \code{ms} is  equal to zero and
#' a general (or biased) random walk otherwise.
#' @seealso \code{\link{sim.Stasis}}, \code{\link{sim.OU}}, \code{\link{as.paleoTS}}
#' @export
#'
#' @examples
#' x.grw <- sim.GRW(ms = 0.5)
#' x.urw <- sim.GRW(ms = 0)
#' plot(x.grw, ylim = range(c(x.grw$mm, x.urw$mm)))
#' plot(x.urw, add = TRUE, col = "blue")
#' legend(x = "topleft", c("GRW", "URW"), col = c("black", "blue"), lty = 1)
sim.GRW <- function (ns = 20, ms = 0, vs=0.1, vp=1, nn=rep(20,ns), tt=0:(ns-1))
  # simulates GRW; ns= number of samples, ms=mean and vs=variance of the step distribution,
  #  vp=population variance, tt=ages of the samples
{
  MM<- array(dim=ns)
  mm<- array(dim=ns)
  vv<- array(dim=ns)
  dt<- diff(tt)

  inc<- rnorm(ns-1, ms*dt, sqrt(vs*dt))	# evolutionary increments

  MM<- cumsum(c(0,inc))	# true means
  mm<- MM + rnorm(ns, 0, sqrt(vp/nn))	# true means plus sampling error
  vv<- rep(vp, ns)

  gp<- c(ms, vs)
  names(gp)<- c("mstep", "vstep")

  res<- as.paleoTS(mm=mm, vv=vv, nn=nn, tt=tt, MM=MM, genpars=gp, label="Created by sim.GRW", reset.time=FALSE)
  return(res)
}


#' Simulate Stasis time-series for trait evolution
#'
#' @param ns number of populations in the sequence
#' @param theta mean of populations
#' @param omega variance among populations
#' @param vp phenotypic variance within populations
#' @param nn vector of population sample sizes
#' @param tt vector of population times (ages)
#'
#' @return a \code{paleoTS} object
#' @seealso \code{\link{sim.GRW}}, \code{\link{sim.OU}}, \code{\link{as.paleoTS}}
#' @export
#'
#' @examples
#' x <- sim.Stasis(omega = 0.5, vp = 0.1)
#' w.sta <- fitSimple(x, model = "Stasis")
#' w.ss  <- fitSimple(x, model = "StrictStasis")
#' compareModels(w.sta, w.ss)
#'
sim.Stasis <- function(ns = 20, theta = 0, omega = 0, vp = 1, nn = rep(20,ns), tt = 0:(ns-1))
  # simulate stasis
{
  xmu<- rnorm(ns, mean=theta, sd=sqrt(omega))
  xobs<- xmu + rnorm(ns, 0, sqrt(vp/nn))
  gp<- c(theta, omega)
  names(gp)<- c("theta", "omega")

  x <- as.paleoTS(mm=xobs,vv=rep(vp,ns),nn=nn,tt=tt,MM=xmu,genpars=gp,label="Created by sim.Stasis", reset.time=FALSE)
  return(x)
}


#' Fit simple models of trait evolution
#'
#' @param y a \code{paleoTS} object
#' @param model the model to be fit, one of \code{"GRW", "URW", "Stasis", "OU",
#'   "covTrack"}
#' @param method parameterization to use: \code{Joint} or \code{AD}; see Details
#' @param pool if TRUE, sample variances are substituted with their pooled
#'   estimate
#' @param z a vector of a covariate, used only for the "covTrack" model
#' @param hess if TRUE, standard errors computed from the Hessian matrix are
#'   returned
#'
#' @return a \code{paleoTSfit} object with the model fitting results
#' @export
#' @importFrom stats dnorm optim rnorm var
#' @details This is a convenience function that calls the specific individual
#' functions for each model and parameterization, such as \code{opt.GRW} and
#' \code{opt.joint.GRW}. The models that this function can fit are:
#' \itemize{
#' \item  \strong{GRW}: General Random Walk. Under this model, evolutionary
#' changes, or "steps" are drawn from a distribution with a mean of \code{mstep}
#' and variance of \code{vstep}.  \code{mstep} determines directionality and
#' \code{vstep} determines volatility (Hunt, 2006).
#' \item  \strong{URW}:
#' Unbiased Random Walk. Same as GRW with \code{mstep} = 0, and thus evolution
#' is non-directional. For a URW, \code{vstep} is the rate parameter.
#' \item \strong{Stasis}: This parameterization follows Sheets & Mitchell (2001), with
#' a constant mean \code{theta} and variance \code{omega} (equivalent to white
#' noise).
#' \item  \strong{Strict Stasis}: Same as Stasis with \code{omega} = 0,
#' indicating no real evolutionary differences; all observed variation is
#' sampling error (Hunt et al. 2015).
#' \item  \strong{OU}: Ornstein-Uhlenbeck
#' model (Hunt et al. 2008). This model is that of a population ascending a
#' nearby peak in the adaptive landscape. The optimal trait value is \code{theta},
#' \code{alpha} indicates the strength of attraction to that peak (= strength of
#' stabilizing selection around \code{theta}), \code{vstep} measures the random walk component (from genetic drift) and \code{anc} is the trait value
#' at the start of the sequence.
#' \item  \strong{covTrack}: Covariate-tracking (Hunt et al. 2010). The trait tracks
#' a covariate with slope \code{b1}, consistent with an adaptive response. \code{evar} is the
#' residual variance, and, under \code{method = "Joint"}, \code{b0} is the intercept of the
#' relationship between trait and covariate.
#' model. }
#' @note For the covariate-tracking model, z should be a vector of length
#' \emph{n} when \code{method = "Joint"} and \emph{n} - 1 when \code{method =
#' "AD"}, where \emph{n} is the number of samples in \code{y}. \cr \cr Method =
#' \code{"Joint"} is a full likelihood approach, considering each time-series as
#' a joint sample from a multivariate normal distribution.  Method = \code{"AD"}
#' is a REML approach that uses the differences between successive samples.
#' They perform similarly, but the Joint approach does better under some
#' circumstances (Hunt, 2008).
#' @references Hunt, G. 2006. Fitting and comparing models of phyletic
#' evolution: random walks and beyond. \emph{Paleobiology} 32(4): 578-601. \cr
#' Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based
#' assessment of modes, rates, punctuations and process. p. 117-131 \emph{In}
#' From Evolution to Geobiology: Research Questions Driving Paleontology at the
#' Start of a New Century. Bambach, R. and P. Kelley (Eds). \cr Hunt, G., M. A.
#' Bell and M. P. Travis. 2008. Evolution toward a new adaptive optimum:
#' phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62(3):
#' 700-710. \cr Sheets, H. D., and C. Mitchell. 2010. Why the null matters:
#' statistical tests, random walks and evolution. \emph{Genetica} 112–
#' 113:105–125. \cr
#'
#' @seealso \code{\link{opt.GRW}}, \code{\link{opt.joint.GRW}},
#'   \code{\link{opt.joint.OU}}, \code{\link{opt.covTrack}}
#' @examples
#' y <- sim.Stasis(ns = 20, omega = 2)
#' w1 <- fitSimple(y, model = "GRW")
#' w2 <- fitSimple(y, model = "URW")
#' w3 <- fitSimple(y, model = "Stasis")
#' compareModels(w1, w2, w3)
fitSimple<- function(y, model = c("GRW", "URW", "Stasis", "StrictStasis", "OU", "covTrack"),
                     method = c("Joint", "AD"), pool = TRUE, z = NULL, hess = FALSE)
{
  model<- match.arg(model)
  method<- match.arg(method)

  if(pool){
    tv<- test.var.het(y)
    pv<- round(tv$p.value, 0)
    wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
    if(pv <= 0.05)	warning(wm)
  }

  if(method=="Joint"){
    if(model=="GRW")			w<- opt.joint.GRW(y, pool=pool, hess=hess)
    if(model=="URW")			w<- opt.joint.URW(y, pool=pool, hess=hess)
    if(model=="Stasis")			w<- opt.joint.Stasis(y, pool=pool, hess=hess)
    if(model=="StrictStasis")	w<- opt.joint.StrictStasis(y, pool=pool, hess=hess)
    if(model=="OU")				w<- opt.joint.OU(y, pool=pool, hess=hess)
    if(model=="covTrack")	{
      if(is.null(z))	stop("Covariate [z] needed for covTrack model.")
      w<- opt.joint.covTrack(y, z, pool=pool, hess=hess)}
  }
  if(method=="AD"){
    if(model=="GRW")			w<- opt.GRW(y, pool=pool, hess=hess)
    if(model=="URW")			w<- opt.URW(y, pool=pool, hess=hess)
    if(model=="Stasis")			w<- opt.Stasis(y, pool=pool, hess=hess)
    if(model=="StrictStasis")	w<- opt.StrictStasis(y, pool=pool, hess=hess)
    if(model=="OU")				stop("AD method not available for OU model.  Consider using Joint method.")
    if(model=="covTrack")	{
      if(is.null(z))	stop("Covariate [z] needed for covTrack model.")
      w<- opt.covTrack(y, z, pool=pool, hess=hess) }
  }

  return(w)
}



#' Fit a set of standard evolutionary models
#'
#' @param y a \code{paleoTS} object
#' @param silent if TRUE, results are returned as a list and not printed
#' @param method "Joint" or "AD", see \code{\link{fitSimple}}
#' @param ... other arguments passed to model fitting functions
#'
#' @details Function \code{fit3models} fits the general (biased) random walk (GRW),
#' unbiased random walk (URW), and Stasis models.  In addition to these three,
#' \code{fit4models} also fits the model of Strict Stasis.
#'
#' @return if silent = FALSE, a table of model fit statistics, also printed to the
#' screen.  if silent = TRUE, a list of the model fit statistics and model parameter values.
#'
#' @seealso \code{\link{fitSimple}}
#' @export
#'
#' @examples
#' x <- sim.GRW(ns = 50, ms = 0.2)
#' fit4models(x)
fit3models<- function (y, silent = FALSE, method = c("Joint", "AD"), ...)
{
  args<- list(...)
  check.var<- TRUE
  if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
  if(all(y$nn ==1)) check.var <- FALSE
  if (check.var){
    tv<- test.var.het(y)
    pv<- round(tv$p.value, 0)
    wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
    if(pv <= 0.05)	warning(wm)
  }

  method <- match.arg(method)
  if (method == "AD") {
    m1 <- opt.GRW(y, ...)
    m2 <- opt.URW(y, ...)
    m3 <- opt.Stasis(y, ...)
  }
  else if (method == "Joint") {
    m1 <- opt.joint.GRW(y, ...)
    m2 <- opt.joint.URW(y, ...)
    m3 <- opt.joint.Stasis(y, ...)
  }
  mc <- compareModels(m1, m2, m3, silent = silent)
  invisible(mc)
}


#' @describeIn fit3models add model of "Strict Stasis" to the three models
#' @export
fit4models<- function(y, silent = FALSE, method = c("Joint", "AD"), ...)
{
  args<- list(...)
  check.var<- TRUE
  if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
  if(all(y$nn ==1)) check.var <- FALSE
  if (check.var){
    tv<- test.var.het(y)
    pv<- round(tv$p.value, 0)
    wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
    if(pv <= 0.05)	warning(wm)
  }

  method <- match.arg(method)
  if (method == "AD") {
    m1 <- opt.GRW(y, ...)
    m2 <- opt.URW(y, ...)
    m3 <- opt.Stasis(y, ...)
    m4 <- opt.StrictStasis(y, ...)
  }
  else if (method == "Joint") {
    m1 <- opt.joint.GRW(y, ...)
    m2 <- opt.joint.URW(y, ...)
    m3 <- opt.joint.Stasis(y, ...)
    m4 <- opt.joint.StrictStasis(y, ...)
  }

  mc <- compareModels(m1, m2, m3, m4, silent = silent)
  invisible(mc)
}


#' Analytical ML estimator for random walk and stasis models
#'
#' @param y a \code{paleoTS} object
#'
#' @return a vector of \code{mstep} and \code{vstep} for \code{mle.GRW},
#' \code{vstep} for \code{mle.URW}, and \code{theta} and \code{omega} for
#' \code{mle.Stasis}
#' @export
#' @note These analytical solutions assume even spacing of samples and equal
#' sampling variance in each, which will usually be violated in real data.
#' They are used here mostly to generate initial parameter estimates for
#' numerical optimization; they not likely to be called directly by the user.
#' @seealso \code{\link{fitSimple}}

mle.GRW<- function(y)
{
  nn<- length(y$mm)-1
  tt<- (y$tt[nn+1]-y$tt[1])/nn
  eps<- 2*pool.var(y)/round(stats::median(y$nn))  # sampling variance
  dy<- diff(y$mm)
  my<- mean(dy)

  mhat<- my/tt
  vhat<- (1/tt)*( (1/nn)*sum(dy^2) - my^2 - eps)

  w<- c(mhat, vhat)
  names(w)<- c("mstep", "vstep")
  return(w)
}

#' @describeIn mle.GRW  ML parameter estimates for URW model
#' @export
mle.URW<- function(y)
  # Gives analytical parameter estimates (URW), assuming:
  #   evenly spaced samples (constant dt)
  #	same sampling variance for each dx (=2*Vp/n)
  # Will try with reasonable values even if assumptions are violated
{
  nn<- length(y$mm)-1
  tt<- (y$tt[nn+1]-y$tt[1])/nn
  eps<- 2*pool.var(y)/round(stats::median(y$nn))  # sampling variance
  dy<- diff(y$mm)
  my<- mean(dy)

  vhat<- (1/tt)*( (1/nn)*sum(dy^2) - eps)

  w<- vhat
  names(w)<- "vstep"
  return(w)
}

#' @describeIn mle.GRW  ML parameter estimates for Stasis model
#' @export
mle.Stasis <- function (y)
  # analytical solution to stasis model
{
  ns<- length(y$mm)
  vp<- pool.var(y)
  th<- mean(y$mm[2:ns])
  om<- var(y$mm[2:ns]) - vp/stats::median(y$nn)

  w<- c(th, om)
  names(w)<- c("theta", "omega")
  return(w)
}


# internal function, not exported
logL.GRW<- function(p, y)
  # function to return log-likelihood of step mean and variance (M,V)= p,
  # given a paleoTS object
{
  # get parameter values: M is Mstep, V is Vstep
  M<- p[1]
  V<- p[2]

  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  S<- dnorm(x=dy, mean=M*dt, sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}

# internal function, not exported
logL.URW<- function(p, y)
  # function to return log-likelihood of step mean and variance (M,V)= p,
  # given a paleoTS object
{
  # get parameter values: V is Vstep
  V<- p[1]

  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  S<- dnorm(x=dy, mean=rep(0,nd), sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}



# internal function, not exported
logL.Stasis <- function(p, y)
  ## logL of stasis model
{
  # get parameter estimates
  M<-p[1]	# M is theta
  V<-p[2]	# V is omega
  dy<- diff(y$mm)
  nd<- length(dy)

  sv<- y$vv/y$nn
  svD<- sv[2:(nd+1)]  # only need sampling variance of descendant
  anc<- y$mm[1:nd]

  S<- dnorm(x=dy, mean=M-anc, sd=sqrt(V + svD), log=TRUE)
  return(sum(S))
}

# internal function, not exported
logL.StrictStasis<- function(p, y)
{
  M <- p[1]
  dy <- diff(y$mm)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svD <- sv[2:(nd + 1)]
  anc <- y$mm[1:nd]
  S <- dnorm(x = dy, mean = M - anc, sd = sqrt(svD), log = TRUE)
  return(sum(S))
}



#' Fit evolutionary model using "AD" parameterization
#'
#' @param y a \code{paleoTS} object
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param meth optimization algorithm, passed to \code{optim()}
#' @param hess if TRUE, return standard errors of parameter estimates from the
#' hessian matrix
#'
#' @return a \code{paleoTSfit} object with the model fitting results
#' @details These functions use differences between consecutive populations in the
#' time series in order to remove temporal autocorrelation.  This is referred to as
#' the "Ancestor-Descendant" or "AD" parameterization by Hunt [2008], and it is a REML
#' approach (like phylogenetic independent contrasts).  A full ML approach, called
#' "Joint" was found to have somewhat better performance (Hunt, 2008) and generally
#' should be used instead.
#' @note  It is easier to use the convenience function \code{fitSimple}.
#' @seealso \code{\link{fitSimple}}, \code{\link{opt.joint.GRW}}
#' @references
#' Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.
#' \emph{Paleobiology} 32(4): 578-601.
#'
#' @export
#'
#' @examples
#' x <- sim.GRW(ns = 20, ms = 1)  # strong trend
#' plot(x)
#' w.grw <- opt.GRW(x)
#' w.urw <- opt.URW(x)
#' compareModels(w.grw, w.urw)
opt.GRW<- function (y, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE)
  # optimize estimates of step mean and variance for GRW model
  # y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.GRW(y)
  if (p0[2] <= 0)	p0[2]<- 1e-7
  names(p0)<- c("mstep", "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)

  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if(inherits(w, "try-error"))
  {
    cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else
      w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if(inherits(w, "try-error"))	# if still fails
    {
      warning("opt.GRW failed ", immediate.=TRUE)
      w$par<- c(NA,NA)
      w$value<- NA
    }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='AD', K=2, n=length(y$mm)-1, se=w$se)

  return (wc)
}

#' @describeIn opt.GRW  fit the URW model by the AD parameterization
#' @export
opt.URW<- function (y, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE)
  # optimize estimates of step mean and variance for GRW model
  # y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.URW(y)
  if (p0 <= 0)	p0<- 1e-7
  names(p0)<- "vstep"
  if (is.null(cl$ndeps))		cl$ndeps<- p0/1e4

  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if(inherits(w, "try-error"))
  {
    cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
    else
      w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if(inherits(w, "try-error"))	# if still fails
    {
      warning("opt.URW failed ", immediate.=TRUE)
      w$par<- NA
      w$value<- NA
    }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='AD', K=1, n=length(y$mm)-1, se=w$se)

  return (wc)
}


#' @describeIn opt.GRW  fit the Stasis model by the AD parameterization
#' @export
opt.Stasis<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
  # optimize estimates of step mean and variance for GRW model
  # y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.Stasis(y)
  if (p0[2] <= 0 || is.na(p0[2]))	p0[2]<- 1e-7
  names(p0)<- c("theta", "omega")
  if (is.null(cl$ndeps))		cl$ndeps<- abs(p0/1e4)

  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if(inherits(w, "try-error"))
  {
    cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else
      w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if(inherits(w, "try-error"))	# if still fails
    {
      warning("opt.Stasis failed ", immediate.=TRUE)
      w$par<- c(NA,NA)
      w$value<- NA
    }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='AD', K=2, n=length(y$mm)-1, se=w$se)

  return (wc)
}


#' @describeIn opt.GRW  fit the Strict Stasis model by the AD parameterization
#' @export
opt.StrictStasis<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
{
  p0 <- mean(y$mm)
  names(p0) <- c("theta")
  if (is.null(cl$ndeps))
    cl$ndeps <- min(abs(p0/10000), 1e-7)
  if (pool)
    y <- pool.var(y, ret.paleoTS = TRUE)
  if (meth == "L-BFGS-B")
    w <- try(optim(p0, fn = logL.StrictStasis, method = meth, lower = c(NA,
                                                                        0), control = cl, hessian = hess, y = y), silent = TRUE)
  else w <- try(optim(p0, fn = logL.StrictStasis, method = meth,
                      control = cl, hessian = hess, y = y), silent = TRUE)
  if(inherits(w, "try-error")) {
    cl$ndeps <- rep(1e-09, length(p0))
    if (meth == "L-BFGS-B")
      w <- try(optim(p0, fn = logL.StrictStasis, method = meth,
                     lower = c(NA, 0), control = cl, hessian = hess,
                     y = y), silent = TRUE)
    else w <- try(optim(p0, fn = logL.StrictStasis, method = meth,
                        control = cl, hessian = hess, y = y), silent = TRUE)
    if(inherits(w, "try-error")) {
      warning("opt.StrictStasis failed ", immediate. = TRUE)
      w$par <- c(NA, NA)
      w$value <- NA
    }
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = "StrictStasis",
                      method = "AD", K = 1, n = length(y$mm) - 1, se = w$se)
  return(wc)
}




# internal function, not exported
logL.joint.GRW<- function (p, y)
  # returns logL of GRW model for paleoTS object x
  # p is vector of parameters: alpha, ms, vs
{
  # prepare calculations
  anc<- p[1]
  ms<- p[2]
  vs<- p[3]
  n<- length(y$mm)

  # compute covariance matrix
  VV<- vs*outer(y$tt, y$tt, FUN=pmin)
  diag(VV)<- diag(VV) + y$vv/y$nn

  # compute logL based on multivariate normal
  M<- rep(anc, n) + ms*y$tt
  S<- mnormt::dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

  return(S)
}

# internal function, not exported
logL.joint.URW<- function(p, y)
  # returns logL of URW model for paleoTS object x
  # p is vector of parameters: anc, vs
{
  # prepare calculations
  anc<- p[1]
  vs<- p[2]
  n<- length(y$mm)

  # compute covariance matrix
  VV<- vs*outer(y$tt, y$tt, FUN=pmin)
  diag(VV)<- diag(VV) + y$vv/y$nn 	# add sampling variance

  # compute logL based on multivariate normal
  M<- rep(anc, n)
  #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
  #S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
  S<- mnormt::dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

  return(S)
}

# internal function, not exported
logL.joint.Stasis<- function (p, y)
  # returns logL of Stasis model for paleoTS object x
  # p is vector of parameters: theta, omega
{
  # prepare calculations
  theta<- p[1]
  omega<- p[2]
  n<- length(y$mm)

  # compute covariance matrix
  VV<- diag(omega + y$vv/y$nn)  # omega + sampling variance

  # compute logL based on multivariate normal
  M<- rep(theta, n)
  #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
  #S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
  S<- mnormt::dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

  return(S)
}

# internal function, not exported
logL.joint.StrictStasis<- function (p, y)
{
  theta<- p[1]
  n<- length(y$mm)
  VV<- diag(y$vv/y$nn)
  detV<- det(VV)
  invV<- solve(VV)
  M<- rep(theta, n)
  #S<- dmvnorm(x$mm, mean = M, sigma = VV, log = TRUE)
  S<- mnormt::dmnorm(y$mm, mean = M, varcov = VV, log = TRUE)

  return(S)
}



#' Fit evolutionary models using the "Joint" parameterization
#'
#' @param y a \code{paleoTS} object
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param meth optimization algorithm, passed to \code{optim()}
#' @param hess if TRUE, return standard errors of parameter estimates from the
#' hessian matrix
#' @return a \code{paleoTSfit} object with the model fitting results
#' @export
#' @details These functions use the joint distribution of population means to fit models
#' using a full maximum-likelihood approach. This approach was found to have somewhat
#' better performance than the "AD" approach, especially for noisy trends (Hunt, 2008).
#' @note  It is easier to use the convenience function \code{fitSimple}.
#' @seealso \code{\link{fitSimple}}, \code{\link{opt.GRW}}
#'
#' @references
#' #' Hunt, G., M. J. Hopkins and S. Lidgard. 2015. Simple versus complex models of trait evolution
#' and stasis as a response to environmental change. PNAS 112(16): 4885-4890.
#'
#' @examples
#' x <- sim.GRW(ns = 20, ms = 1)  # strong trend
#' plot(x)
#' w.grw <- opt.joint.GRW(x)
#' w.urw <- opt.joint.URW(x)
#' compareModels(w.grw, w.urw)
opt.joint.GRW<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
  # optimize GRW model using alternate formulation
{
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  ## get initial estimates
  p0<- array(dim=3)
  p0[1]<- y$mm[1]
  p0[2:3]<- mle.GRW(y)
  if (p0[3]<=0)	p0[3]<- 1e-7
  names(p0)<- c("anc", "mstep", "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  cl$ndeps[cl$ndeps==0]<- 1e-8

  if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, lower=c(NA,NA,0), hessian=hess, y=y)
  else 					w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, hessian=hess, y=y)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='Joint', K=3, n=length(y$mm), se=w$se)

  return (wc)
}

#' @describeIn opt.joint.GRW  fit the URW model by the Joint parameterization
#' @export
opt.joint.URW<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
  # optimize URW model using alternate formulation
{
  ## check if pooled
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  ## get initial estimates
  p0<- array(dim=2)
  p0[1]<- y$mm[1]
  p0[2]<- min(c(mle.URW(y), 1e-7)) ## handles negative vstep estimates
  names(p0)<- c("anc","vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  cl$ndeps[cl$ndeps==0]<- 1e-8

  if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, lower=c(NA,0), hessian=hess, y=y)
  else					w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, hessian=hess, y=y)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='Joint', K=2, n=length(y$mm), se=w$se)

  return (wc)
}

#' @describeIn opt.joint.GRW  fit the Stasis model by the Joint parameterization
#' @export
opt.joint.Stasis<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
  # optimize Stasis model using alternate formulation
{
  ## check if pooled
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  ## get initial estimates
  p0<- mle.Stasis(y)
  if(p0[2]<=0 || is.na(p0[2]))	p0[2]<- 1e-7
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  cl$ndeps[cl$ndeps==0]<- 1e-9

  if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, lower=c(NA,0), hessian=hess, y=y)
  else				  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, hessian=hess, y=y)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='Joint', K=2, n=length(y$mm), se=w$se)

  return (wc)

}


#' @describeIn opt.joint.GRW  fit the Strict Stasis model by the Joint parameterization
#' @export
opt.joint.StrictStasis<- function (y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE)
{
  if (pool)
    y <- pool.var(y, ret.paleoTS = TRUE)
  p0 <- mean(y$mm)
  names(p0)<- "theta"
  w <- optim(p0, fn = logL.joint.StrictStasis, control = cl, method = "Brent", lower=min(y$mm), upper=max(y$mm),
             hessian = hess, y = y)
  names(w$par)<- "theta"
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="StrictStasis", method="Joint", K=1, n=length(y$mm), se=w$se)
  return(wc)
}


# internal functions, not exported
ou.M<- function(anc, theta, aa, tt) theta*(1 - exp(-aa*tt)) + anc*exp(-aa*tt)
ou.V<- function(vs, aa, tt)        (vs/(2*aa))*(1 - exp(-2*aa*tt))


#' Simulate an Ornstein-Uhlenbeck time-series
#'
#' @param ns  number of populations in the sequence
#' @param anc ancestral phenotype
#' @param theta OU optimum (long-term mean)
#' @param alpha strength of attraction to the optimum
#' @param vstep step variance
#' @param vp phenotypic variance of each sample
#' @param nn vector of sample sizes
#' @param tt vector of sample times (ages)
#'
#' @return a \code{paleoTS} object
#' @details This function simulates an Ornstein-Uhlenbeck (OU) process. In
#'   microevolutionary terms, this models a population ascending a nearby peak
#'   in the adaptive landscape. The optimal trait value is \code{theta},
#'   \code{alpha} indicates the strength of attraction to that peak (= strength
#'   of stabilizing selection around \code{theta}), \code{vstep} measures the
#'   random walk component (from genetic drift) and \code{anc} is the trait
#'   value at the start of the sequence.
#' @references Hunt, G., M. A. Bell and M. P. Travis. 2008. Evolution toward a new adaptive
#' optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62(3):
#' 700-710.
#' @examples
#' x1 <- sim.OU(alpha = 0.8)  # strong alpha
#' x2 <- sim.OU(alpha = 0.1)  # wearker alpha
#' plot(x1)
#' plot(x2, add = TRUE, col = "blue")
#'
#' @export
#'
#' @seealso \code{\link{opt.joint.OU}}
sim.OU<- function (ns = 20, anc = 0, theta = 10, alpha = 0.3, vstep = 0.1, vp = 1, nn = rep(20, ns),
                   tt = 0:(ns-1))
  ## generate a paleoTS sequence according to an OU model
{

  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  MM[1] <- anc
  x <- rnorm(nn[1], mean = MM[1], sd = sqrt(vp))
  mm[1] <- mean(x)
  vv[1] <- var(x)
  for (i in 2:ns) {
    ex<- ou.M(MM[i-1], theta, alpha, dt[i-1])
    vx<- ou.V(vstep, alpha, dt[i-1])
    MM[i]<- rnorm(1, ex, sqrt(vx))
    x <- rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
    mm[i] <- mean(x)
    vv[i] <- var(x)
  }

  gp <- c(anc, theta, alpha, vstep)
  names(gp) <- c("anc", "theta", "alpha", "vstep")
  res <- as.paleoTS(mm = as.vector(mm), vv = as.vector(vv), nn = nn, tt = tt, MM = MM,
                    genpars = gp, label = "Created by sim.OU()", reset.time=FALSE)

  return(res)
}

# internal function, not exported
logL.joint.OU<- function(p, y)
  # returns logL of OU model for paleoTS object x
  #
{
  # prepare calculations
  anc<- p[1]
  vs<- p[2]
  theta<- p[3]
  aa<- p[4]
  n<- length(y$mm)

  # compute covariance matrix
  ff<- function (a,b) abs(a-b)
  VV<- outer(y$tt, y$tt, FUN=ff)
  VV<- exp(-aa*VV)
  VVd<- ou.V(vs,aa,y$tt)
  VV2<- outer(VVd,VVd,pmin)
  VV<- VV*VV2
  diag(VV)<- VVd + y$vv/y$nn 	# add sampling variance
  #detV<- det(VV)
  #invV<- solve(VV)

  # compute logL based on multivariate normal
  M<- ou.M(anc, theta, aa, y$tt)
  #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
  #S<- dmvnorm(t(x$mm), mean=M, sigma=VV, log=TRUE)
  S<- mnormt::dmnorm(t(y$mm), mean=M, varcov=VV, log=TRUE)
  return(S)
}


#' Fit Ornstein-Uhlenbeck model using the "Joint" parameterization
#'
#' @param y a \code{paleoTS} object
#' @param pool if TRUE, sample variances are substituted with their pooled estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param meth optimization algorithm, passed to \code{optim()}
#' @param hess if TRUE, return standard errors of parameter estimates from the
#' hessian matrix
#' @return a \code{paleoTSfit} object with the model fitting results
#' @export
#' @details This function fits an Ornstein-Uhlenbeck (OU) model to time-series data. The OU
#' model has four generating parameters: an ancestral trait value (\code{anc}), an optimum
#' value (\code{theta}), the strength of attraction to the optimum (\code{alpha}), and a
#' parameter that reflects the tendency of traits to diffuse (\code{vstep}).  In a
#' microevolutionary context, these parameters can be related to natural selection and
#' genetic drift; see Hunt et al. (2008).
#' @note  It is easier to use the convenience function \code{fitSimple}.  Note also that
#' preliminary work found that the "AD" parameterization did not perform as well for the OU
#' model and thus it is not implemented here.
#' @seealso \code{\link{fitSimple}}, \code{\link{opt.joint.GRW}}
#' @references Hunt, G., M. A. Bell and M. P. Travis. 2008. Evolution toward a new adaptive
#' optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62(3):
#' 700-710.
#'
#' @examples
#' x <- sim.OU(vs = 0.5)  # most defaults OK
#' w <- opt.joint.OU(x)
#' plot(x, modelFit = w)
opt.joint.OU<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
{
  ## check if pooled
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  ## get initial estimates
  w0<- mle.GRW(y)
  halft<- (y$tt[length(y$tt)]-y$tt[1])/4			# set half life to 1/4 of length of sequence
  p0<- c(y$mm[1], w0[2]/10, y$mm[length(y$mm)], log(2)/halft)
  names(p0)<- c("anc","vstep","theta","alpha")
  #print(p0)
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  cl$ndeps[cl$ndeps==0]<- 1e-8


  if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, lower=c(NA,1e-10,NA,1e-8), hessian=hess, y=y)
  else 				  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, hessian=hess, y=y)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='OU', method='Joint', K=4, n=length(y$mm), se=w$se)

  return (wc)
}




#' Simulate trait evolution that tracks a covariate
#'
#' @param ns number of populations in a sequence
#' @param b slope of the relationship between the change in the covariate and
#'   the change in the trait
#' @param evar residual variance of the same relationship
#' @param z vector of covariate that the trait tracks
#' @param nn vector of sample sizes for populations
#' @param tt vector of times (ages) for populations
#' @param vp phenotypic trait variance within each population
#'
#' @details In this model, changes in a trait are linearly related to changes in
#'   a covariate with a slope of \code{b} and residual variance \code{evar}:
#'   \code{dx = b * dz + eps}, where \code{eps ~ N(0, evar)}.  This model was
#'   described, and applied to an example in which body size changes tracked
#'   changes in temperature, by Hunt et al. (2010).
#'
#' @note  For a trait sequence of length \code{ns}, the covariate, \code{z}, can
#'   be of length \code{ns} - 1,in which case it is interpreted as the vector of
#'   \emph{changes}, \code{dz}. If \code{z} is of length \code{ns},
#'   differences are taken and these are used as the \code{dz}'s.
#'
#' @return a \code{paleoTS} object
#' @export
#' @references Hunt, G, S. Wicaksono, J. E. Brown, and K. G. Macleod. 2010.
#' Climate-driven body size trends in the ostracod fauna of the deep Indian
#' Ocean. \emph{Palaeontology} 53(6): 1255-1268.
#' @examples
#' set.seed(13)
#' z <- c(1, 2, 2, 4, 0, 8, 2, 3, 1, 9, 4, 3)
#' x <- sim.covTrack(ns = 12, z = z, b = 0.5, evar = 0.2)
#' plot(x, ylim = c(-1, 10))
#' lines(x$tt, z, col = "blue")
sim.covTrack<- function (ns=20, b=1, evar=0.1, z, nn=rep(20, times=ns), tt=0:(ns-1), vp=1)
{
  # check on length of covariate
  if(length(z)==ns)
  {
    z<- diff(z)
    warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
  }

  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- rnorm(ns-1, b*z, sqrt(evar))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))  # add sampling error
  vv <- rep(vp, ns)
  gp <- c(b, evar)
  names(gp) <- c("slope.b", "evar")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM,
                    genpars = gp, label = "Created by sim.covTrack()")

  return(res)
}



# internal function, not exported
logL.covTrack<- function(p, y, z)
  # z is covariate; pars = b (slope), evar (variance)
  # IMPT: z is of length ns-1; one for each AD transition IN ORDER
{
  b<- p[1]
  evar<- p[2]
  dy <- diff(y$mm)
  dt <- diff(y$tt)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svA <- sv[1:nd]
  svD <- sv[2:(nd + 1)]
  svAD <- svA + svD
  #S <- -0.5 * log(2 * pi * (evar + svAD)) - ((dy - (b*z))^2)/(2 * (evar + svAD))
  S<- dnorm(x=dy, mean=b*z, sd=sqrt(evar+svAD), log=TRUE)
  return(sum(S))
}





#' Fit a model in which a trait tracks a covariate
#'
#' @param y a \code{paloeTS} object
#' @param z a vector of covariate values
#' @param pool if TRUE, sample variances are substituted with their pooled
#'   estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param meth optimization algorithm, passed to \code{optim()}
#' @param hess if TRUE, return standard errors of parameter estimates from the
#'   hessian matrix

#' @details In this model, changes in a trait are linearly related to changes in
#'   a covariate with a slope of \code{b} and residual variance \code{evar}:
#'   \code{dx = b * dz + eps}, where \code{eps ~ N(0, evar)}.  This model was
#'   described, and applied to an example in which body size changes tracked
#'   changes in temperature, by Hunt et al. (2010). \cr
#'
#'   For the AD version (\code{opt.covTrack}), a trait sequence of
#'   length \code{ns}, the covariate, \code{z}, can be of length \code{ns} - 1,
#'   interpreted as the vector of \emph{changes}, \code{dx}. If \code{z} is
#'   of length \code{ns}, differences are taken and these are used as the
#'   \code{dx}'s, with a warning issued. \cr
#'
#'   The Joint version
#'   (\code{opt.joint.covTrack}), \code{z} should be of length \code{ns} and
#'   there is an additional parameter that is the intercept of the linear
#'   relationship between trait and covariate. See warning below about using the
#'   Joint version.
#'
#' @section Warning: The Joint parameterization of this model can be fooled by
#'   temporal autocorrelation and, especially, trends in the trait and the
#'   covariate.  The latter is tested for, but the AD parameterization is
#'   generally safer for this model.

#' @return a \code{paleoTSfit} object with the results of the model fitting
#' @export

#' @references Hunt, G, S. Wicaksono, J. E. Brown, and K. G. Macleod. 2010. Climate-driven
#' body size trends in the ostracod fauna of the deep Indian Ocean. \emph{Palaeontology}
#' 53(6): 1255-1268.
#' @seealso \code{\link{fitSimple}}
#' @examples
#' set.seed(13)
#' z <- c(1, 2, 2, 4, 0, 8, 2, 3, 1, 9, 4, 3)
#' x <- sim.covTrack(ns = 12, z = z, b = 0.5, evar = 0.2)
#' w.urw <- opt.URW(x)
#' w.cov <- opt.covTrack(x, z = z)
#' compareModels(w.urw, w.cov)
opt.covTrack<- function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
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
  reg<- stats::lm(diff(y$mm) ~ z-1)
  p0<- c(stats::coef(reg), var(stats::resid(reg)))
  names(p0) <- c("b", "evar")

  # pool variances if needed and do optimization
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (is.null(cl$ndeps))
    cl$ndeps <- abs(p0/10000)
  cl$ndeps[cl$ndeps==0]<- 1e-8  ## will fail o.w. if any p0=0
  if (meth == "L-BFGS-B")
    w <- optim(p0, fn=logL.covTrack, method = meth, lower = c(NA, 0), control = cl, hessian = hess, y=y, z=z)
  else  w<- optim(p0, fn=logL.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)


  if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  return(wc)
}



# internal function, not exported
logL.joint.covTrack<- function(p, y, z)
  # z is covariate; pars = b0 (intercept), b1 (slope), evar (variance)
  # z is of length ns; one for each sample in the time-series
{
  b0<- p[1]
  b1<- p[2]
  evar<- p[3]
  sv <- y$vv/y$nn
  S<- dnorm(x=y$mm, mean=b0 + b1*z, sd=sqrt(evar+sv), log=TRUE)
  return(sum(S))
}

#' @describeIn opt.covTrack  fits the covTrack model using the joint parameterization
#' @export
opt.joint.covTrack<- function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
{
  # check if z is of proper length; first difference if necessary
  ns<- length(y$mm)
  if (length(z) != ns)  stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n" )

  # check if both covariate and trait are trended; if so, warn that this approach is probably unreliable
  r.trait<- stats::cor(y$mm, y$tt)
  r.cov<- stats::cor(z, y$tt)
  mess<- c("Both the trait and covariate are strongly trended.  The joint approach is probably unreliable in this situation;
           consider using the AD parameterization instead.")
  if(abs(r.trait) > 0.6 && abs(r.cov > 0.6))	warning(mess)


  # get initial estimates by regression
  reg<- stats::lm(y$mm ~ z)
  p0<- c(stats::coef(reg), var(stats::resid(reg)))
  names(p0) <- c("b0", "b1", "evar")

  # pool variances if needed and do optimization
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (is.null(cl$ndeps))
    cl$ndeps <- rep(1e-8, length(p0))
  if (meth == "L-BFGS-B")
    w <- optim(p0, fn=logL.joint.covTrack, method = meth, lower = c(NA,NA,0), control = cl, hessian = hess, y=y, z=z)
  else  w<- optim(p0, fn=logL.joint.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)


  if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='Joint', K=3, n=length(y$mm), se=w$se)
  return(wc)
}


