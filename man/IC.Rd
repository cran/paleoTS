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
  \item{n}{ sample size for AICc and BIC calculations (= the number of evolutionary transisitons) }
  \item{meth}{ which information criterion to compute; one of AIC, AICc, or BIC }
  \item{aa}{ vector of AIC or AICc values used to compute Akaike weights }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
 the computed information criterion, or a vector of Akaike weights
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}} } 

\examples{
 x <- sim.GRW(ns=40, ms=0.1, vs=0.1)
 m1<- opt.GRW(x)	
 m2<- opt.URW(x)
 m3<- opt.Stasis(x)
 akaike.wts(c(m1$AICc, m2$AICc, m3$AICc))  # note easier to use fit3models()
	
}

\keyword{models}
\keyword{ts}