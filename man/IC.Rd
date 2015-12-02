\name{IC}
\alias{IC}
\alias{akaike.wts}

\title{ Compute information criterion scores and Akaike weights for evoltuionary models }
\description{
 These functions compute information criteria (\code{IC}) or Akaike weights based on information scores (\code{akaike.wts}).
 \cr
 Function \code{IC} is used internally and generally will not need to be called directly by the user.
}
\usage{
IC(logL, K, n = NULL, method = c("AICc", "AIC", "BIC"))
akaike.wts(aa)
}

\arguments{
  \item{logL}{ log-likelihood }
  \item{K}{ the number of free parameters }
  \item{n}{ sample size for AICc and BIC calculations }
  \item{method}{ which information criterion to compute; one of AIC, AICc, or BIC }
  \item{aa}{ vector of AIC or AICc values used to compute Akaike weights }
}
\details{
These functions are used by the functions \code{as.paleoTSfit} and \code{compareModels}, and will rarely need to be used directly by the user.}
  
\value{
 the computed information criterion, or a vector of Akaike weights
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}}, \code{\link{compareModels}}, \code{\link{as.paleoTSfit}} } 

\examples{
 x <- sim.GRW(ns=40, ms=0.1, vs=0.1)
 m1<- opt.GRW(x)	
 m2<- opt.URW(x)
 aw<- akaike.wts(c(m1$AICc, m2$AIC))  # easier to use compareModels(m1,m2)
	
}

\keyword{models}
\keyword{ts}