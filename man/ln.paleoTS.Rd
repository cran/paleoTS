% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basics.R
\name{ln.paleoTS}
\alias{ln.paleoTS}
\title{Approximate log-transformation of time-series data}
\usage{
ln.paleoTS(y)
}
\arguments{
\item{y}{a \code{paleoTS} object}
}
\value{
the converted \code{paleoTS} object
}
\description{
Approximate log-transformation of time-series data
}
\details{
For a random variable x, its approximate mean on a natural log
  scale is the log of its untransformed mean.  The approximate variance on a
  log scale is equal to the squared coefficient of variation.
}
\note{
This transformation only makes sense for variables with dimension and a
  true zero point, such as lengths and areas.
}
\examples{
x <- sim.Stasis(ns = 10, theta = 20, omega = 1)
plot(x)
xl <- ln.paleoTS(x)
plot(xl)
}
\references{
Hunt, G. 2006. Fitting and comparing models of phyletic
  evolution: random walks and beyond.  \emph{Paleobiology} 32:578-601. \cr
  Lewontin, R. 1966. On the measurement of relative variability.
  \emph{Systematic Zoology} 15:141-142.
}
