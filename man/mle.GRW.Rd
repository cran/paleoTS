\name{mle.GRW}
\alias{mle.GRW}
\alias{mle.URW}
\alias{mle.Stasis}

\title{ Maximum likelihood parameter estimators }
\description{
  These functions calculate maximum likelihood estimators for the general random walk (\code{mle.GRW}), unbiased random walk (\code{mle.URW}) and stasis (\code{mle.Stasis}) models.
}

\usage{
mle.GRW(y)
mle.URW(y)
mle.Stasis(y)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
}
\details{
  These functions return maximum likelihood estimators for the general random walk \code{mle.GRW}, unbiased random walk \code{mle.URW} and stasis \code{mle.Stasis} models, but \strong{only under a restriced set of circumstances are these valid!}  For these estimators to be valid, the sampling error must be the same in all samples, which generally means equal sample size and variances in all samples.  For the random walk models, it is also assumed that samples are evenly spaced in time.  Because these assumptions usually do not hold for paleontological data, almost all users should instead use the numerical optimization functions (see \code{\link{opt.GRW}}).  The main purpose for the present functions is to provide starting estimates for numerical optimization. 

}
\value{
  a vector of parameter estimates, either \code{c(mstep, vstep)}, or \code{c(theta, omega)}
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}}, \code{\link{logL.GRW}}  }

\examples{
  ## Warning: better to use opt.RW() for real data
  y <- sim.GRW(ns=20, ms=0, vs=1)
  w <- mle.GRW(y)
}
\keyword{models}
\keyword{ts}