\name{mle.rw}
\alias{mle.rw}
\alias{mle.stasis}

\title{ Maximum likelihood parameter estimators }
\description{
  These functions calculate maximum likelihood estimators for the general random walk (\code{mle.rw}) and stasis (\code{mle.stasis}) models.
}

\usage{
mle.rw(x)
mle.stasis(x)
}

\arguments{
  \item{x}{ a \code{paleoTS} object }
}
\details{
  These functions return maximum likelihood estimators for the general random walk \code{mle.rw} and stasis \code{mle.stasis} models, but \strong{only under a restriced set of circumstances are these valid!}.  For these estimators to be valid, the sampling error must be the same in all samples, which generally means equal sample size and variances in all samples.  For the general random walk model, it is also assumed that samples are evenly spaced in time.  Because these assumptions usually do not hold for paleontological data, almost all users should instead use the numerical optimization functions (see \code{\link{opt.RW}}).  The main purpose for the present functions is to provide starting estimates for numerical optimization. 

}
\value{
  a vector of parameter estimates, either \code{c(mstep, vstep)}, or \code{c(theta, omega)}
}

\references{ Hunt, G. 2007. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}, in press. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.RW}}, \code{\link{logL.RW}}  }

\examples{
  ## Warning: better to use opt.RW() for real data
  y <- sim.rw(ns=20, 0, 1)
  w <- mle.rw(y)
}
\keyword{models}
\keyword{ts}