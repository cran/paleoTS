\name{paleoTS-package}
\alias{paleoTS-package}
\alias{paleoTS}
\docType{package}
\title{
Analysis of evoltuionary time-series 
}
\description{
This package includes a variety of functions useful for analyzing evoltuionary time-series, 
particularly for paleontological data.  The functions are useful for analyzing sequences in terms 
of statistical models, with three models implemented: directional evolution (modeled as a general random walk),
unbiased random walks, and stasis.  These analyses are outlined in Hunt (2006).
}
\details{
\tabular{ll}{
Package: \tab paleoTS\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2007-10-10\cr
License: \tab GPL 2\cr
}

}
\author{ Gene Hunt}

\references{
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601.
}

\keyword{ package }
\examples{
 x <- sim.GRW(ns=30, ms=0, vs=.3)	# simulate unbiased random walk
 fit3models(x)		# compare fits of directional, random walk, and stasis models
}
