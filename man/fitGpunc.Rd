\name{fitGpunc}
\alias{fitGpunc}
\alias{opt.punc}
\alias{opt.joint.punc}
\alias{logL.punc}
\alias{logL.punc.omega}
\alias{logL.joint.punc}
\alias{logL.joint.punc.omega}


\title{ Analyze evolutionary models with unsampled punctuations }

\description{
 Functions required to fit evolutionary models with puntuations that are rapid relative to the temporal spacing 
 of samples (so-called unsampled punctuations).
}

\usage{
fitGpunc(y, ng = 2, minb = 7, pool = TRUE, oshare = TRUE, method=c('Joint', 'AD'), 
         silent = FALSE, hess=FALSE, parallel = FALSE, ...)
         
opt.punc(y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", 
         hess = FALSE, oshare)
opt.joint.punc(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", 
               hess=FALSE, oshare)
logL.punc(p, y, gg)
logL.punc.omega(p, y, gg)
logL.joint.punc(p, y, gg)
logL.joint.punc.omega(p, y, gg)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{ng}{ the number of separate segments in the sequence }
  \item{minb}{ the minimum number of samples within a segment to consider }
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{oshare}{ logical, if \code{TRUE}, the same variance (\code{omega}) is assumed across all segments.  If \code{FALSE}, separate variances are assumd for each segment  }
  \item{method}{ parameterization to use: based on ancestor-descendant (AD) differences, or on Joint consideration of all samples}
  \item{silent}{ if \code{TRUE}, less information is printed to the screen as the model is fit }
  \item{hess}{ if \code{TRUE}, standard errors are computed from the Hessian matrix}
  \item{parallel}{ logical, if TRUE, the analysis is parallelized using \code{foreach}}
  \item{\dots}{ other arguments to send to \code{opt.punc}}
  \item{p}{ parameters of the punctuation model to be optimized}
  \item{gg}{ numeric vector indicating membership of each sample in segments 1, 2, .. ng}
  \item{cl}{ control list to be passed to \code{optim}}
  \item{meth}{ optimization method, to be passed to \code{optim} } 
}
\details{
  These functions are used to fit a model with an unsampled punctuation.
  It is equivalent to a Stasis model in which the optimum instantaneously shifts at one or more points in time; see references below for details.
  Users are likely only to use \code{fitGpunc}, which will calls the other functions in order to find the best parameter
  estimates and shift points for the segments.  
  
  The \code{parallel} option uses \code{doParallel} and \code{foreach}.  It detects the number of available cores and uses that number minus one for the parallelization.
}
\value{
 The log-likelihood functions return the log-likelihood of the model for a given set of parameter values (\code{p}), 
 assuming that the periods of Stasis have the same variance (\code{logL.punc.omega}) or different variances (\code{logL.punc}).
 
 Functions \code{fitGpunc} and \code{opt.punc} return a \code{paleoTSfit} object.

In addition, the \code{paleoTSfit} object returned by function \code{fitGpunc} has the following additional elements:
  \item{all.logl}{ log-likelihoods for all tested partitions of the series into segments}
  \item{GG}{ matrix of indices of initial samples of each tested segment configuration; each column of \code{GG} corresponds to the elements of \code{all.logl}}
}

\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601.  \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology} \bold{ 34}:360--377.  \cr
}

\author{ Gene Hunt }

\seealso{ \code{\link{sim.punc}}, \code{\link{opt.GRW}}, \code{\link{fit.sgs}} }

\examples{
x<- sim.punc(theta=c(0,5), ns=c(20,20), omega=c(0.5,0.5), vp=0.2)
w<- fitGpunc(x, ng=2, minb=7, pool=TRUE, oshare=TRUE)
print (w$parameters)
## plot using modelFit argument to show the solution
plot(x, modelFit=w)

## now, plot logL with respect to shift point
shift.times<- x$tt[w$GG]
plot(shift.times, w$all.logl, typ="o", cex=0.5, 
     main="Log-likelihood of different shift times")

}

\keyword{models}
\keyword{ts}