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
fitGpunc(y, ng = 2, minb = 5, pool = TRUE, oshare = TRUE, method=c('AD', 'Joint'), silent = FALSE, hess=FALSE, ...)
opt.punc(y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE, oshare)
opt.joint.punc(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare)
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
}
\value{
 The log-likelihood functions return the log-likelihood of the model for a given set of parameter values (\code{p}), 
 assuming that the periods of Stasis have the same variance (\code{logL.punc.omega}) or different variances (\code{logL.punc}).
 
 Functions \code{fitGpunc} and \code{opt.punc} return a list with the following elements:
  \item{par }{ parameter estimates}
  \item{value }{ the log-likelihood of the optimal solution}
  \item{counts}{ returned by \code{optim} }
  \item{convergence}{ returned by \code{optim} }
  \item{message}{ returned by \code{optim} }
  \item{p0}{ initial guess for parameter values at start of optimization}
  \item{K}{ number of parameters in the model}
  \item{n}{ the number of observations, equal to the number of evoltuionary transistions}
  \item{AIC}{ Akaike information criterion}
  \item{AICc}{ modified Akaike information criterion}
  \item{BIC}{ Bayes (or Schwarz) information criterion}
  \item{se}{ standard errors for parameter estimates, computed from the curvature of the log-likelihood surface (only if \code{hess = TRUE})}
  \item{\dots}{ other output from call to \code{optim}}

In addition, function \code{fitGpunc} also returns the following elements:
  \item{shift.start}{ index of the last sample before a new segment}
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
x<- sim.punc(theta=c(0,5), ns=c(20,20), omega=c(0.5,0.5), vp=c(0.2,0.2))
w<- fitGpunc(x, ng=2, minb=7, pool=TRUE, oshare=TRUE)
print (w$parameters)
## plot using modelFit argument to show the solution
plot(x, modelFit=w)
}

\keyword{models}
\keyword{ts}