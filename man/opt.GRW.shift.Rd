\name{opt.GRW.shift}
\alias{opt.GRW.shift}

\title{ Functions for random walks with shifting parameters }

\description{
Functions to simulate and to infer a model with random walk dynamics, with parameter values that 
shift at one or more points in the sequence.
}

\usage{
opt.GRW.shift(y, ng = 2, minb = 7, model = 1, pool = TRUE, silent = FALSE)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{ng}{ the number of different segments in the sequence }
  \item{minb}{ the minimum number of samples to consider as a segment }
  \item{model}{ options for variants of random walk to fit (see \code{Details}). }
  \item{pool}{ if TRUE, pool phenotypic variances across samples }
  \item{silent}{ if TRUE, do not print information on fitting to screen }
}
\details{
This model divides an evolutionary sequence into two or more non-overlapping parts 
called segments, and then fits some version of the random walk model to each segment separately. \cr
The \code{model} argument has four options: \cr
  \code{model=1}, general random walk model, step variance shared across segments \cr
  \code{model=2}, general random walk model, step mean shared across segments \cr
  \code{model=3}, unbiased random walk, separate step variance for each segment
  \code{model=4}, general random walk, speparate step mean and variance for each segment
}
\value{
  A \code{paleoTSfit} object.
}
\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601.  \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology} \bold{ 34}:360--377.  \cr
}
\author{ Gene Hunt}

\seealso{ \code{\link{sim.GRW.shift}}, \code{\link{opt.GRW}}, \code{\link{fitMult}}, \code{\link{as.paleoTSfit}} }
\examples{
 x<- sim.GRW.shift(ns=c(20,20), ms=c(0,1), vs=c(0.2, 0.2))
 plot(x)
 m.urw<- fitSimple(x, model="URW", method="AD")
 m.shift<- opt.GRW.shift(x, model=1)	# only AD method is available for this model
 compareModels(m.urw, m.shift)
 shift.time<- x$tt[m.shift$parameters["shift1"]]
 abline(v=shift.time, col="red", abline=3)
}

\keyword{models}
\keyword{ts}