\name{opt.RW}
\alias{opt.RW}
\alias{opt.RWu}
\alias{opt.stasis}

\title{ Numerically find maximum likelihood solutions }

\description{
  Functions to find maximum likelihood solutions to general random walk (\code{opt.RW}), unbiased random walk \code{opt.RWu}, and stasis models \code{opt.stasis}.
}

\usage{
opt.RW(y, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
opt.RWu(y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.stasis(y, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{cl}{ control list, passed to function \code{optim} }
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{meth}{ optimization method, passed to function \code{optim} }
  \item{hess}{ logical, indicating whether to calculate standard errors from the Hessian matrix }
}
\details{
  These functions numerically search a log-likelihood surface for its optimum--they are a convenient wrapper to \code{optim}.  Arguments \code{meth}, \code{cl}, and \code{hess} are passed to \code{optim}; see that function's help for details.  These are included to allow sophisticated users greater control over the optimization; most users will just use the defaults. For \code{meth}, only \code{meth="L-BFGS-B"} and \code{meth="BFGS"} are implemented.  The former constrains estimates of the step variance to be non-negative, which can be useful paramters which cannot truly be negative, such as vstep (random walk) and omega (stasis model). \cr

  Initial estimates to start the optimization come from analytical solutions based on assuming equal sampling error across samples and evenly spaced samples in time (functions \code{mle.rw} and \code{mle.stasis}). 
  
}

\value{
  A list including:
  \item{par }{ parameter estimates}
  \item{value }{ the log-likelihood of the optimal solution}
  \item{counts}{ returned by \code{optim} }
  \item{convergence}{ returned by \code{optim} }
  \item{message}{ returned by \code{optim} }
  \item{p0}{ initial guess for parameter values at start of optimization}
  \item{K}{ number of parameters in the model}
  \item{AIC}{ Akaike information criterion}
  \item{AICc}{ modified Akaike information criterion}
  \item{BIC}{ Bayes (or Schwarz) information criterion}
}

\references{ Hunt, G. 2007. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\note{
  Standard errors computed from the Hessian matrix are reasonably accurate for \code{mstep} and \code{theta}, but not useful for the vstep and omega because of the asymmetry of the log-likelihood surfaces. 
}

\seealso{\code{\link{logL.RW}}, \code{\link{fit3models}}, \code{\link{opt.RW.mult}}, \code{\link{sim.rw}} }

\examples{
 ## generate data for a directional sequence
 y <- sim.rw(30, sm=1, sv=1)
 plot(y)
 m.rw<- opt.RW(y)
 m.rwu<- opt.RWu(y)
 m.sta<- opt.stasis(y)

 ## print log-likelihoods, easier to use function fit3models()
 cat(m.rw$value, m.rwu$value, m.sta$value, "\n")
}

\keyword{models}
\keyword{ts}