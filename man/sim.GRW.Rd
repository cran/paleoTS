% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpleModels.R
\name{sim.GRW}
\alias{sim.GRW}
\title{Simulate random walk or directional time-series for trait evolution}
\usage{
sim.GRW(ns = 20, ms = 0, vs = 0.1, vp = 1, nn = rep(20, ns), tt = 0:(ns - 1))
}
\arguments{
\item{ns}{number of populations in the sequence}

\item{ms}{mean of evolutionary "steps"}

\item{vs}{variance of evolutionary "steps"}

\item{vp}{phenotypic variance within populations}

\item{nn}{vector of population sample sizes}

\item{tt}{vector of population times (ages)}
}
\value{
a \code{paleoTS} object
}
\description{
Simulate random walk or directional time-series for trait evolution
}
\details{
The general random walk model considers time in discrete steps.
At each time step, an evolutionary change is drawn at random from a distribution of
possible evolutionary "steps."  It turns out that the long-term dynamics of an evolving
lineage depend only on the mean and variance of this step distribution.  The former,
\code{mstep}, determined the directionality in a sequence and the latter, \code{vstep},
determines its volatility.
}
\note{
This function simulates an unbiased random walk if \code{ms} is  equal to zero and
a general (or biased) random walk otherwise.
}
\examples{
x.grw <- sim.GRW(ms = 0.5)
x.urw <- sim.GRW(ms = 0)
plot(x.grw, ylim = range(c(x.grw$mm, x.urw$mm)))
plot(x.urw, add = TRUE, col = "blue")
legend(x = "topleft", c("GRW", "URW"), col = c("black", "blue"), lty = 1)
}
\seealso{
\code{\link{sim.Stasis}}, \code{\link{sim.OU}}, \code{\link{as.paleoTS}}
}
