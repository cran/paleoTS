\name{fit3models}
\alias{fit3models}

\title{Do model fits for three evolutionary models}

\description{
  Function fits three models to an evolutionary time series: (1) general random walk, (2) unbiased random walk, and (3) stasis. 
}

\usage{
fit3models(y, pool = TRUE, silent = FALSE, wts = "AICc")
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{pool}{ logical, if TRUE, variances are pooled across samples }
  \item{silent}{ logical, if TRUE, results are not printed }
  \item{wts}{ what version of AIC to use for calculation of Akaike weights; either \code{AIC} or \code{AICc} }
}

\details{
  
}
\value{
  If \code{silent=FALSE}, function \code{fit3models} returns a vector of parameter estimates; if \code{silent=TRUE}, a list with the following components is returned, with models listed in the order of general random walk, unbiased random walk, and stasis: 
  \item{aic }{Akaike information criterion}
  \item{aicc }{AIC modified for low sample sizes}
  \item{logl}{log-likelihoods of the three models}
  \item{hats}{parameter estimates for general random walk (mstep, vstep), unbiased random walk (vstep) and stasis (theta, omega) models}
  \item{ak.wts}{vector of Akaike weights}
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}} }

\examples{
 y<- sim.GRW(ns=30, ms=0.5, vs=0.2)
 fit3models(y)
 y<- sim.Stasis(ns=30, theta=0, omega=0)
 fit3models(y)
}
\keyword{models}
\keyword{ts}