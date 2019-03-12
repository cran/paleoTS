# rates functions #

#'@title Log-rate, Log-interval (LRI) method of Gingerich
#'@description Gingerich (1993) introduced a method that plots on log-log scale,
#'  the rate and interval for each pair of samples in an evolutionary sequence.
#'  On this plot, the slope is interpreted as an indicator of evolutionary mode
#'  (-1 for stasis, 0.5 for random walk, 0 for directional), and the intercept
#'  is interpreted as a measure of the rate of evolution over one generation.
#'
#'@param y a \code{paleoTS} object
#'@param gen.per.t the number of generations per unit time
#'@param draw logical, if TRUE, a plot is produced
#'
#'@return A named vector with three elements: \code{Intercept}, \code{slope}, and
#'  \code{GenerationalRate}
#'@export
#'@details Following Gingerich (1993), a robust line is fit through the
#'  points by minimizing the sum of absolute deviations. If generations are one
#'  year long and time is measured in Myr, \code{gen.per.t}= 1e6.
#'@note This method was important in early attempts to disentangle
#'  evolutionary tempo and mode. I view likelihood-based methods as more
#'  informative, and in particular the estimation of 'Generational Rates' using
#'  LRI is compromised by sampling error; see Hunt (2012) and the example below.
#'@references
#'  Gingerich, P.D. 1993. Quantification and comparison of
#'  evolutionary rates. \emph{American Journal of Science} 293-A:453–478. \cr \cr
#'  Hunt, G. 2012. Measuring rates of phenotypic evolution and the
#'  inseparability of tempo and mode. \emph{Paleobiology} 38:351–373.
#'@seealso \code{\link{lynchD}}
#' @examples
#' set.seed(1)
#' xFast <- sim.GRW(ns = 20, ms = 0.5, vs = 0.2)  # fast evolution
#' xSlow <- sim.Stasis(ns = 20, omega = 0)        # strict stasis (zero rates)
#' lri.Fast <- LRI(xFast, draw = FALSE)
#' lri.Slow <- LRI(xSlow, draw = FALSE)
#' print(lri.Fast[3], 4)
#' print(lri.Slow[3], 4)  # LRI thinks strict stasis rates are faster!
LRI<- function(y, gen.per.t=1e6, draw=TRUE)
{
  n<- length(y$mm)
  y$tt<- y$tt*gen.per.t  # convert to generational timescale (e.g., 1e6 when time is measured in Myr and generation time = 1 yr)
  sp<- sqrt(pool.var(y))
  dy<- outer(y$mm, y$mm, FUN='-')  # outer() much faster than looping
  dt<- outer(y$tt, y$tt, FUN='-')
  dy<- abs(dy)
  dt<- abs(dt)

  dy<- dy/sp
  ut<- upper.tri(dy)
  logI<- log10(dt[ut])
  logR<- log10(dy[ut]) - logI

  # subset only non-zero rates
  ok<- is.finite(logR)
  if (sum(ok) != length(logI))	warning("Some zero rates were ignored.")

  # call robust regression function
  lf<- opt.lad(logI[ok], logR[ok])
  lf[3]<- 10^lf[1]   # generational rate
  names(lf)<- c("Intercept", "slope", "GenerationalRate")

  # do LRI plot, if desired
  if (draw)
  {
    xl<- c(0, max(logI[ok]))
    yl<- c(min(logR[ok]), max(lf[1], max(logR[ok])))
    plot(logI[ok], logR[ok], xlim=xl, ylim=yl, xlab='log10 Interval [generations]', ylab='log10 Rate', cex=0.6)
    abline(lf[1], lf[2], col='black', lwd=2)
    title("LRI plot")
    mtext(paste('data label: ', y$lab), cex=0.6, font=3, col='darkgrey')
    restext<- paste('Slope = ', round(lf[2],3), '\n', 'Intercept = ', round(lf[1],3), '\n', 'Generational Rate = ', round(lf[3],5), '\n', sep='')
    text(0, min(logR[ok]), restext, adj=c(0,0), cex=0.7, font=2)
  }

  #w<- list(b0=b0, b1=b1, h0=h0, logR=logR, logI=logI, dy=dy)
  return(lf)
}

# internal function, not exported
# function to fit line in LRI plot robustly, using min abs deviations as per Gingerich 1993
opt.lad<- function (x, y)
{
  ok<- is.finite(x) & is.finite(y)
  xok<- x[ok]
  yok<- y[ok]
  w.ls<- stats::lm(yok ~ xok)  # use LS coef as starting point in optimization

  # function minimized for least abs deviation
  fad<- function(p, x, y)  sum( abs(y - p[2]*x - p[1]) )
  w.lad<- stats::optim(w.ls$coef, fn=fad, x=xok, y=yok)
  return(w.lad$par)
}



#' @title Compute Lynch's Delta rate metric
#' @description This function computes D, the rate metric proposed by Lynch
#'   (1990).  This metric derives from the random walk model, with \code{D =
#'   Vstep/(2Vp)}, where \code{Vstep} is the step variance of the unbiased
#'   random walk, and \code{Vp} is the within sample variance, pooled among
#'   samples.  Under mutation - drift equilibrium, \code{D} is expected to range
#'   approximately between 5e-5 and 5e-3.
#'
#' @param y a \code{paleoTS} object
#' @param gen.per.t the number of generations per unit time
#' @param pool logical, whether variances should be pooled over samples
#' @param method parameterization to use: based on ancestor-descendant (AD)
#'   differences, or Joint consideration of all samples
#' @param ... further arguments, passed to \code{opt.URW} or \code{opt.joint.URW}
#'
#' @return
#' \item{D}{  value of rate metric }
#' \item{pooled.var}{ value
#'   of pooled within-sample variance }
#'   \item{gen.per.t}{ number of generations per unit time }
#'   \item{vstep}{ computed \code{Vstep}, at the original time
#'   scale of \code{y}}
#'   \item{drift.range}{ expected minimum and maximum values
#'   of \code{D} consistent with neutral evolution }
#'   \item{result}{ conclusion reached about the plausibility of neutral evolution }
#' @export
#' @references
#'  Lynch (1990). The rate of morphological evolution in mammals from the
#'  standpoint of the neutral expectation. \emph{The American Naturalist}
#'  136:727-741.
#'   Hunt, G. 2012. Fitting and comparing models of phyletic evolution: random
#'   walks and beyond.  \emph{Paleobiology} 38:351-373.
#' @examples
#' y <- sim.GRW(ns = 20, ms = 0, vs = 1e-4, tt=seq(0, 1e6, length.out=20))  # per-year simulation
#' lynchD(y, gen.per.t = 1)
lynchD<- function (y, gen.per.t=1e6, pool=TRUE, method=c('Joint', 'AD'), ...)
{
  method<- match.arg(method)
  vp<- pool.var(y)
  if (method == 'AD')
  {
    wu<- opt.URW(y, pool=pool)
    vs<- unname(wu$par)
  }
  if (method == 'Joint')
  {
    wu<- opt.joint.URW(y, pool=pool)
    vs<- unname(wu$par[2])
  }

  D<- 0.5*(vs/vp)/gen.per.t

  drift.min<- 5e-5
  drift.max<- 5e-3

  if (D < drift.min)	res<- "Slower than drift range"
  else if (D > drift.max)	res<- "Faster than drift range"
  else					res<- "Within range of drift"

  w<- list(D=D, pooled.var=vp, gen.per.t=gen.per.t, vstep=vs, drift.range=c(drift.min, drift.max), result=res)
  return(w)
}






#' @title Compute Expected Squared Divergence (ESD) for Evolutionary Models
#' @description  Computes for a specified model and duration of time the
#'   expected squared divergence (ESD), which is a useful measure of the
#'   magnitude or rate of change across different models.
#'
#' @param y a \code{paleoTS} object
#' @param dt the time interval to evaluate ESD
#' @param model the model of evolution to assume; see Details
#' @param method Joint or AD parameterization
#' @param pool logical, if TRUE, variances are averaged (pooled) across samples
#' @param ... other arguments to the model-fitting functions
#'
#' @return the ESD value
#' @export
#'
#' @details Hunt (2012) argued that rate metrics make sense only in the context
#'   of specific models of evolution.  It is thus difficult to meaningfully
#'   compare rates across sequences generated by different evolutionary
#'   processes.  ESD values can be used for a specified model and duration as a
#'   comparable measure of the amount of evolutionary change that is expected.
#'   Acceptable values for the model argument can be "GRW" for the general
#'   random walk (directional change), "URW" for the unbiased random walk, and
#'   "Stasis."  In addition, one can also specify "allThree", in which case all
#'   these models will be fit and the resulting ESD will be the weighted average
#'   of them, using model support (Akaike weights) for the weighting (see Hunt
#'   [2012], p. 370)
#' @references
#' Hunt, G. 2012. Measuring rates of phenotypic evolution and the
#'  inseparability of tempo and mode. \emph{Paleobiology} 38:351–373.
#' @examples
#'   x<- sim.GRW(ns=20)
#'   esd.urw<- ESD(x, dt=10, model="URW")
#'   esd.all<- ESD(x, dt=10, model="allThree")
ESD<- function(y, dt, model=c("GRW", "URW", "Stasis", "allThree"), method=c("Joint", "AD"), pool=TRUE, ...)
{
  model<- match.arg(model)
  method<- match.arg(method)

  if(model=="GRW"){
    if(method=="Joint")		w<- opt.joint.GRW(y, pool=pool, ...)
    else if (method=="AD")	w<- opt.GRW(y, pool=pool, ...)
    pp<- w$par
    esd<- pp["mstep"]^2 * dt*2 + pp["vstep"] * dt
  }
  if(model=="URW"){
    if(method=="Joint")		w<- opt.joint.URW(y, pool=pool, ...)
    else if (method=="AD")	w<- opt.URW(y, pool=pool, ...)
    pp<- w$par
    esd<- pp["vstep"] * dt

  }
  if(model=="Stasis"){
    if(method=="Joint")		w<- opt.joint.Stasis(y, pool=pool, ...)
    else if (method=="AD")	w<- opt.Stasis(y, pool=pool, ...)
    pp<- w$par
    esd<- 2*pp["omega"]

  }
  if(model=="allThree"){
    w<- fit3models(y, silent=TRUE, method=method, pool=pool)
    aw<- w$modelFits$Akaike.wt
    mn<- rownames(w$modelFits)

    esd.grw<- ESD(y, dt, model="GRW", method=method, pool=pool)
    esd.urw<- ESD(y, dt, model="URW", method=method, pool=pool)
    esd.sta<- ESD(y, dt, model="Stasis", method=method, pool=pool)

    esdV<- c(esd.grw, esd.urw, esd.sta)
    esd<- sum(esdV * aw)
  }

  return(unname(esd))

}
