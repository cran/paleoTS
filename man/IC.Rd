\name{IC}
\alias{IC}
\alias{akaike.wts}

\title{ Compute information criterion scores and Akaike weights for evoltuionary models }
\description{
 These functions compute information criteria (\code{IC}) or Akaike weights based on information scores (\code{akaike.wts}).
 \cr
 Function \code{IC} is used internally by the optimization functions and generally will not need to be called directly by the user.
}
\usage{
IC(w, logL = NULL, K = NULL, n = NULL, meth = c("AICc", "AIC", "BIC"))
akaike.wts(aa)
}

\arguments{
  \item{w}{ output from an optimization such as \code{opt.GRW} }
  \item{logL}{ log-likelihood }
  \item{K}{ the number of free parameters }
  \item{n}{ sample size for AICc and BIC calculations }
  \item{meth}{ which information criterion to compute; one of AIC, AICc, or BIC }
  \item{aa}{ vector of AIC or AICc values used to compute Akaike weights }
}
\details{
Function \code{IC} can take two kinds of arguments.  The easiest is to send it the output of any of the \code{paleoTS} 
optimizations, in which case the function will automaically extract the necessary information.  Alternatively, the 
log-likelihoods, number of parameters and sample size can be passed explicitly.  \cr
Function \code{akaike.wts} takes a vector of AIC or AICc values and computes corresponding Akaike weights.
}
  
\value{
 the computed information criterion, or a vector of Akaike weights
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}} } 

\examples{
 x <- sim.GRW(ns=40, ms=0.1, vs=0.1)
 m1<- opt.GRW(x)	
 m2<- opt.URW(x)
 m3<- opt.Stasis(x)
 akaike.wts(c(m1$AICc, m2$AICc, m3$AICc))  # it is easier to use fit3models()
	
}

\keyword{models}
\keyword{ts}