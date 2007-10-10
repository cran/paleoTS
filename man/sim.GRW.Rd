\name{sim.GRW}
\alias{sim.GRW}
\alias{sim.Stasis}

\title{ Simulate evolutionary time-series }
\description{
  Simulate the evolution of a trait according to general random walk or stasis models.
}

\usage{
sim.GRW(ns = 20, ms = 0, vs = 0.1, vp = 1, nn = rep(20, ns), tt = 1:ns)
sim.Stasis(ns = 20, theta = 0, omega = 0, vp = 1, nn = rep(20,ns), tt = 1:ns)
}

\arguments{
  \item{ns}{ number of samples in time-series }
  \item{ms}{ mean of the step distribution, random walk model }
  \item{vs}{ variance of the step distribution, random walk model }
  \item{vp}{ within-population trait variance }
  \item{nn}{ vector of the number of individuals in each sample }
  \item{tt}{ vector of sample ages, increases from oldest to youngest }
  \item{theta}{ evolutionary optimum, stasis model }
  \item{omega}{ evolutionary variance, stasis model }
}
\details{
  See reference below for details on parameterization of the models.  \cr
  Briefly, the general random walk model considers time in discrete steps.  The duration of steps does not matter as long as many steps occur between sampled populations.  At each time step, an evolutionary change is drawn at random from a distribution of possible evolutionary "steps."  It turns out that the long-term dynamics of an evolving lineage depend only on the mean and variance of this step distribution.  The former, \code{mstep}, determined the directionality in a sequence and the latter, \code{vstep}, determines its volatility.  \cr
  The stasis model is based on the parameterization of Sheets and Mitchell (2001).  Under this model, there is an evolutionary optimum, \code{theta}, with some amount of true variance, \code{omega}, around this optimum.  
}
\value{
  A \code{paleoTS} object.
}
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. 
\cr
Sheets, H. D., and C. E. Mitchell. 2001. Why the null matters: statistical tests, random walks and evolution. \emph{Genetica} \bold{112-113}:105-125.
}
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}} }
\examples{
 ## generate and plot two paleoTS objects
 y.rw <- sim.GRW(ns=20, ms=0.5, vs=0.1)  
 y.st <- sim.Stasis(ns=20)
 layout(1:2)
 plot(y.rw, col="red")
 plot(y.st, col="blue")
 layout(1)
} 

\keyword{models}
\keyword{ts}