\name{fitMult}

\alias{fitMult}
\alias{opt.joint.Mult}
\alias{opt.Mult}
\alias{opt.RW.SameMs}
\alias{opt.RW.SameVs}
\alias{logL.joint.Mult}
\alias{logL.Mult}
\alias{logL.SameMs}
\alias{logL.SameVs}

\title{ Functions to estimate models over multiple time-series }

\description{
  These functions are used to estimate parameters of models over two or more sequences.  Models that can be fit include the general (directional) random walk, unbiased random walk and stasis, and variants of the general random walk
  in which the step variance but not the mean step is shared over sequences (\code{opt.RW.SameVs}) or the mean step but not the step variance is shared (\code{opt.RW.SameMs}).
}
\usage{
fitMult(yl, model=c("GRW", "URW", "Stasis", "covTrack"), method=c("Joint", "AD"), 
        pool=TRUE, zl=NULL, hess=FALSE)

opt.joint.Mult(yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), 
               pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.Mult(yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), pool=TRUE, 
         meth="L-BFGS-B", hess=FALSE)
opt.RW.SameMs(yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.RW.SameVs(yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
logL.joint.Mult(p, yl, model=c("GRW", "URW", "Stasis"))
logL.Mult(p, yl, model = c("GRW", "URW", "Stasis"))
logL.SameMs(p, yl)
logL.SameVs(p, yl)
}

\arguments{
  \item{yl}{ a list of \code{paleoTS} objects. }
  \item{model}{ \code{GRW} for general random walk (directional evolution), \code{URW} for unbiased random walk, or \code{Stasis} }
  \item{method}{ parameterization to use: \code{Joint} or \code{AD};  see Details}
  \item{pool}{ logical, if TRUE variances are pooled across samples }
  \item{zl}{ a list of covariate vectors, one for each \code{paleoTS} object in \code{yl} }
  \item{cl}{ optimization option, passed to \code{optim} }
  \item{meth}{ optimization option, passed to \code{optim} }
  \item{hess}{ optimization option, passed to \code{optim} }  
  \item{p}{ a vector of parameter values }
}
\details{
  Users will generally only call \code{fitMult}, which can fit general random walks, unbiased random walks or stasis.  Functions \code{opt.RW.SameVs} and \code{opt.RW.SameMs} fit a variant of the general random walk in which only the mean step (or the step variance) is shared across sequences; see Hunt (2006, p. 590) and note that only the \code{AD} method is available for these.
  
  These functions work just as their conterparts for the analysis of single sequences; see those help functions for more detail.
}

\value{
Varies by function, see corresponding functions for the analysis of single sequences for more information.
}

\references{ 
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. \cr
	Hunt, et al. 2010. Climate-driven body-size trends in the ostracod fauna of the deep Indian Ocean, \emph{Palaeontology} \bold{ 53}:1255--1268. 
}

\author{ Gene Hunt }

\seealso{  \code{\link{logL.GRW}}, \code{\link{opt.GRW}}, \code{\link{fitSimple}}, \code{\link{opt.covTrack}} }

\examples{
 ## create two sequences, with different parameter values
 y1<- sim.GRW(ns=20, ms=0, vs=1)
 y2<- sim.GRW(ns=20, ms=0, vs=0.2)

 ## fit some models with at least some shared dynamics across sequences
 m1<- fitMult(list(y1,y2), model="GRW", method="Joint")
 m2<- fitMult(list(y1,y2), model="URW", method="Joint")
 m3<- fitMult(list(y1,y2), model="Stasis", method="Joint")
 compareModels(m1, m2, m3)
 
}

\keyword{models}
\keyword{ts}