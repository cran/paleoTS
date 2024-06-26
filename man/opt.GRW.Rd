% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpleModels.R
\name{opt.GRW}
\alias{opt.GRW}
\alias{opt.URW}
\alias{opt.Stasis}
\alias{opt.StrictStasis}
\title{Fit evolutionary model using "AD" parameterization}
\usage{
opt.GRW(
  y,
  pool = TRUE,
  cl = list(fnscale = -1),
  meth = "L-BFGS-B",
  hess = FALSE
)

opt.URW(
  y,
  pool = TRUE,
  cl = list(fnscale = -1),
  meth = "L-BFGS-B",
  hess = FALSE
)

opt.Stasis(
  y,
  pool = TRUE,
  cl = list(fnscale = -1),
  meth = "L-BFGS-B",
  hess = FALSE
)

opt.StrictStasis(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE)
}
\arguments{
\item{y}{a \code{paleoTS} object}

\item{pool}{if TRUE, sample variances are substituted with their pooled estimate}

\item{cl}{optional control list, passed to \code{optim()}}

\item{meth}{optimization algorithm, passed to \code{optim()}}

\item{hess}{if TRUE, return standard errors of parameter estimates from the
hessian matrix}
}
\value{
a \code{paleoTSfit} object with the model fitting results
}
\description{
Fit evolutionary model using "AD" parameterization
}
\details{
These functions use differences between consecutive populations in the
time series in order to remove temporal autocorrelation.  This is referred to as
the "Ancestor-Descendant" or "AD" parameterization by Hunt [2008], and it is a REML
approach (like phylogenetic independent contrasts).  A full ML approach, called
"Joint" was found to have generally better performance (Hunt, 2008) and generally
should be used instead.
}
\section{Functions}{
\itemize{
\item \code{opt.URW()}: fit the URW model by the AD parameterization

\item \code{opt.Stasis()}: fit the Stasis model by the AD parameterization

\item \code{opt.StrictStasis()}: fit the Strict Stasis model by the AD parameterization

}}
\note{
It is easier to use the convenience function \code{fitSimple}.
}
\examples{
x <- sim.GRW(ns = 20, ms = 1)  # strong trend
plot(x)
w.grw <- opt.GRW(x)
w.urw <- opt.URW(x)
compareModels(w.grw, w.urw)
}
\references{
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.
\emph{Paleobiology} 32(4): 578-601.
}
\seealso{
\code{\link{fitSimple}}, \code{\link{opt.joint.GRW}}
}
