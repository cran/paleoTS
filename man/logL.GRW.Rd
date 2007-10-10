\name{logL.GRW}
\alias{logL.GRW}
\alias{logL.URW}
\alias{logL.Stasis}

\title{ Compute log-likelihoods for random walk and stasis models }

\description{
  Returns log-likelihood for general random walk \code{logL.GRW}, unbiased random walk \code{logL.URW}, and stasis \code{logL.Stasis} models.
}
\usage{
logL.GRW(p, y)
logL.URW(p, y)
logL.Stasis(p, y)
}

\arguments{
  \item{p}{ vector of parameters }
  \item{y}{ a \code{paleoTS} object }
 }
\details{
  For the general random walk, \code{p = c(mstep, vstep)}; for an unbiased random walk, \code{p = vstep}; for the stasis model, \code{p = c(theta, omega)}.  In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values. 
}

\value{
  The log-likelihood of the parameter estimates (\code{p}), given the data (\code{y}).
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }


\seealso{ \code{\link{mle.GRW}}, \code{\link{opt.GRW}} }

\examples{
y<- sim.GRW(20, 0, 1)
L1 <- logL.GRW(p=c(0,1), y)   # actual parameters
L2 <- logL.GRW(p=c(10,10), y) # should be a bad guess
cat (L1, L2, "\n")
}

\keyword{models}
\keyword{ts}