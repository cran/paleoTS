\name{logL.RW}
\alias{logL.RW}
\alias{logL.RWu}
\alias{logL.stasis}

\title{ Compute log-likelihoods for random walk and stasis models }

\description{
  Returns log-likelihood for general random walk \code{logL.RW}, unbiased random walk \code{logL.RWu}, and stasis \code{logL.stasis} models.
}
\usage{
logL.RW(p, y, pool = TRUE)
logL.RWu(p, y, pool = TRUE)
logL.stasis(p, y, pool = TRUE)
}

\arguments{
  \item{p}{ vector of parameters }
  \item{y}{ a \code{paleoTS} object }
  \item{pool}{ logical, if TRUE, variances are pooled across samples }
}
\details{
  For the general random walk, \code{p = c(mstep, vstep)}; for an unbiased random walk, \code{p = vstep}; for the stasis model, \code{p = c(theta, omega)}.  In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values. 
}

\value{
  The log-likelihood of the parameter estimates (\code{p}), given the data (\code{y}).
}

\references{ Hunt, G. 2007. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}, in press. }
\author{ Gene Hunt }


\seealso{ \code{\link{mle.rw}}, \code{\link{opt.RW}} }

\examples{
y<- sim.rw(20, 0, 1)
L1 <- logL.RW(p=c(0,1), y)   # actual parameters
L2 <- logL.RW(p=c(10,10), y) # should be a bad guess
cat (L1, L2, "\n")
}

\keyword{models}
\keyword{ts}