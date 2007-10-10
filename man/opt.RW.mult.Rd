\name{opt.RW.Mult}
\alias{opt.RW.Mult}
\alias{logL.Mult}
\alias{logL.SameMs}
\alias{logL.SameVs}
\alias{opt.RW.SameMs}
\alias{opt.RW.SameVs}

\title{ Functions to analyze multiple time-series jointly }

\description{
  These functions are used to estimate parameters of the general random walk over two or more sequences.  Three different models may be fit: the same parameters are present in all time series (\code{opt.RW.mult}); the same mean step but different step variances in each series (\code{opt.RW.SameMs}); the same step variance but different mean steps in each series (\code{opt.RW.SameVs}).
}
\usage{
logL.Mult(p, y, model = c("GRW", "URW"))
logL.SameMs(p, y)
logL.SameVs(p, y)
opt.RW.Mult(y, cl=list(fnscale=-1), model=c("GRW", "URW"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.RW.SameMs(y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.RW.SameVs(y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
}

\arguments{
  \item{p}{ a vector of parameter values }
  \item{y}{ a emph{list} of \code{paleoTS} objects. }
  \item{model}{ \code{RW} for general random walk, \code{RWu} for unbiased random walk }
  \item{pool}{ logical, if TRUE variances are pooled across samples }
  \item{cl}{ optimization option, passed to \code{optim} }
  \item{meth}{ optimization option, passed to \code{optim} }
  \item{hess}{ optimization option, passed to \code{optim} }
}
\details{
  These functions work just as their conterparts for the analysis of single sequences; see those help functions for more detail.
}

\value{
Varies by function, see corresponding functions for the analysis of single sequences for more information.
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{  \code{\link{logL.GRW}}, \code{\link{opt.GRW}}, \code{\link{fit3models}} }

\examples{
 ## create two sequences, with different parameter values
 y1<- sim.GRW(ns=20, ms=0, vs=1)
 y2<- sim.GRW(ns=20, ms=0, vs=0.2)

 ## fit some models with at least some shared dynamics across sequences
 m1<- opt.RW.Mult(list(y1,y2), model="GRW")
 m2<- opt.RW.SameMs(list(y1, y2))
 m3<- opt.RW.SameVs(list(y1,y2))

 ## fit separate models to each sequence
 msep1<- opt.GRW(y1)
 msep2<- opt.GRW(y2)
 sep.logL<- msep1$value + msep2$value   # total logL sums across two solutions
 sep.aic <- -2*sep.logL + 2*4     # AIC of separate solutions

 ## check out resulting log-likelihoods and AIC values
 mres<- array(dim=c(4,2))
 row.names(mres)<- c("m1", "m2", "m3", "msep")
 colnames(mres)<- c("logL", "AIC")
 mres[1,]<- c(m1$value, m1$AIC)
 mres[2,]<- c(m2$value, m2$AIC)
 mres[3,]<- c(m3$value, m3$AIC)
 mres[4,]<- c(sep.logL, sep.aic)
 print (mres)
}

\keyword{models}
\keyword{ts}