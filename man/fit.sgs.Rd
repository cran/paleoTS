\name{fit.sgs}
\alias{fit.sgs}
\alias{opt.sgs}
\alias{logL.sgs}
\alias{logL.sgs.omega}

\title{ Analyze evolutionary models with well-sampled punctuations }

\description{
 Functions required to fit evolutionary models with sampled puntuations, i.e., where the transitional period 
 is represented by at least several sampled populations.  
}

\usage{
fit.sgs(y, minb = 5, oshare = TRUE, pool = TRUE, silent = FALSE, hess = FALSE, meth = "L-BFGS-B", model = "GRW")
opt.sgs(y, gg, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE, oshare = TRUE, model = "GRW")
logL.sgs(p, y, gg, model = "GRW")
logL.sgs.omega(p, y, gg, model = "GRW")
}

\arguments{
  \item{y}{ a \code{paleoTS} object}
  \item{minb}{ the minimum number of samples within a segment to consider}
  \item{oshare}{ logical, if \code{TRUE}, the same variance (\code{omega}) is assumed across the starting and ending Stasis segments.  If \code{FALSE}, separate variances are assumed}
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{silent}{ if \code{TRUE}, less information is printed to the screen as the model is fit }
  \item{hess}{ if \code{TRUE}, standard errors are computed from the Hessian matrix }
  \item{meth}{ optimization method, to be passed to \code{optim} }
  \item{model}{ either \code{GRW} or \code{URW}, indicating whether evolution during the transitional interval is directional (general random walk) or not (unbiased random walk) }
  \item{p}{ parameters of the punctuation model for the log-likelihood functions }
  \item{gg}{ numeric vector indicating membership of each sample in segments 1, 2, .. ng}
  \item{cl}{ control list to be passed to \code{optim}}
}

\details{
  These functions are used to fit a model with an sampled punctuation.
  Formally, this is a three-segment model that starts as Stasis, transitions to a punctuation of directional evolution (general random walk) or
  unconstrained (unbiased random walk). The name comes from an abbreviation of the three modes in the segments: \emph{S}tasis - \emph{G}eneral Random Walk - \emph{S}tasis,
  bearing in mind that the general random walk can be changed to an unbiased random walk.
  Users are likely only to use \code{fit.sgs}, which will calls the other functions in order to find the best parameter
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

In addition, function \code{fit.sgs} also returns the following elements:
  \item{shift.start}{ index of each sample that starts a new segment}
  \item{all.logl}{ log-likelihoods for all tested partitions of the series into segments}
  \item{GG}{ matrix of indices of initial samples of each tested segment configuration; each column of \code{GG} corresponds to the elements of \code{all.logl}}
}

\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601.  \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology}\bold{34}:In press.  \cr
	}
	
\author{ Gene Hunt }


\seealso{ \code{\link{sim.sgs}}, \code{\link{opt.GRW}}, \code{\link{fitGpunc}}  }

\examples{
 x<- sim.sgs(ns=c(15, 10, 15), ms=0.5, vs=0.3)
 plot(x)
 # compare sampled punctuation to uniform models
 w1<- fit.sgs(x, minb=7, model="GRW")
 wu<- fit3models(x, silent=TRUE)
 aa<- akaike.wts(c(w1$AICc, wu$aicc))
 names(aa)[1]<- "Samp.Punc"
 cat("Akaike Weights:\n")
 print(round(aa,5))
}

\keyword{models}
\keyword{ts}