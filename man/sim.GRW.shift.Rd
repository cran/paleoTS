% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexModels.R
\name{sim.GRW.shift}
\alias{sim.GRW.shift}
\title{Simulate (general) random walk with shift(s) in generating parameters}
\usage{
sim.GRW.shift(
  ns = c(10, 10),
  ms = c(0, 1),
  vs = c(0.5, 0.5),
  nn = rep(30, sum(ns)),
  tt = 0:(sum(ns) - 1),
  vp = 1
)
}
\arguments{
\item{ns}{vector of the number of samples in each segment}

\item{ms}{vector of mean step parameter in each segment}

\item{vs}{vector of step variance parameter in each segment}

\item{nn}{vector of sample sizes, one for each population}

\item{tt}{vector of samples times (ages)}

\item{vp}{phenotypic variance in each sample}
}
\value{
a \code{paleoTS} object with the simulated time-series
}
\description{
Simulate (general) random walk with shift(s) in generating parameters
}
\details{
Simulates under a model in which a sequence is divided into two or more segments.
Trait evolution proceeds as a general random walk, with each segment getting its own
generating parameters (\code{mstep}, \code{vstep}).
}
\examples{
x <- sim.GRW.shift(ns = c(10,10,10), ms = c(0, 1, 0), vs = c(0.1,0.1,0.1))
plot(x)
abline(v = c(9.5, 19.5), lty = 3, lwd = 2, col = "blue")  # shows where dynamics shift
text (c(5, 15, 25), c(2,2,2), paste("segement", 1:3, sep =" "), col = "blue")
}
\seealso{
\code{\link{sim.GRW}}, \code{\link{sim.sgs}}, \code{\link{opt.GRW.shift}}
}
