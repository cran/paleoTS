\name{bootSimpleComplex}
\alias{bootSimpleComplex}

\title{
Use parametric bootstrapping to test the fit of a complex model relative to a simpler one
}
\description{
Tests if a complex model fit is significantly better than that of a simple model, treating the simple model as the null hypothesis and using the log-likelihood ratio as a test statistic.
}
\usage{
bootSimpleComplex(y, simpleFit, complexFit, nboot = 99, minb = 7, 
                  ret.full.distribution = FALSE, parallel = FALSE, ...)
}

\arguments{
  \item{y}{ a \code{paleoTS} object}
  \item{simpleFit}{ a \code{paleoTSfit} object corresponding to a simple model (GRW, URW, Stasis or StrictStasis)}
  \item{complexFit}{ a \code{paleoTSfit} object corresponding to a complex model (from fitGpunc or fitStasisRW)}
  \item{nboot}{ number of boostrap replications}
  \item{minb}{ the minimum number of segments to require before and after the mode shift}
  \item{ret.full.distribution}{ logical, if TRUE, the vectore of test statistics for all replicates is returned}
  \item{parallel}{ logical, if TRUE, the replicates are parallelized using \code{foreach}}
  \item{\dots}{ further arguments, passed to the optimization functions}
}
\details{
Simulations performed in Hunt et al. (2015) showed that the AICc can unduly favor complex models in this context.  The behavior is known in other situations in which a large number of possible shifts are evaluated.

The \code{parallel} option uses \code{doParallel} and \code{foreach}.  It detects the number of available cores and uses that number minus one for the parallelization.
}
\value{
Returns a \code{list}:
\item{LRobs}{ the observed log-likelihood ratio statistic}
\item{p.value}{ p-value of the test}
\item{nullLR}{ optional, if \code{ret.full.distribution}=TRUE, the vector of LR test statistics for the bootstrap replications}
}
\references{
Hunt, G., M. J. Hopkins, and S. L. Lidgard 2015. Simple versus complex models of trait evolution and stasis as a response to environmental change.  \emph{PNAS} \bold{ 112}:4885--4890. \cr
}
\author{
Gene Hunt
}



\seealso{
\code{\link{fit9models}}, \code{\link{fitModeShift}}
}
\examples{
	x<- sim.GRW(ns=40)
	mS<- fitSimple(x, model="URW", method="Joint")
	mC<- fitModeShift(x, order="Stasis-RW", rw.model="URW", method="Joint")
	compareModels(mS, mC)  ## AICc comparisons
	
	\donttest{
	boot<- bootSimpleComplex(x, mS, mC, nboot=20)
	cat("Bootstrap p-value: ", boot$p.value, "\n")
	}
}
\keyword{ models }
\keyword{ ts }
