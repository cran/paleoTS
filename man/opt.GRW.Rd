\name{opt.GRW}
\alias{opt.GRW}
\alias{opt.URW}
\alias{opt.Stasis}

\title{ Numerically find maximum likelihood solutions to evolutionary models }

\description{
  Functions to find maximum likelihood solutions to general random walk (\code{opt.GRW}), unbiased random walk \code{opt.URW}, and stasis models \code{opt.Stasis}.
}

\usage{
opt.GRW(y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.URW(y, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE)
opt.Stasis(y, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE)
}


\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{cl}{ control list, passed to function \code{optim} }
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{meth}{ optimization method, passed to function \code{optim} }
  \item{hess}{ logical, indicating whether to calculate standard errors from the Hessian matrix }
}
\details{
  These functions numerically search a log-likelihood surface for its optimum--they are a convenient wrapper to \code{optim}.  Arguments \code{meth}, \code{cl}, and \code{hess} are passed to \code{optim}; see that function's help for details.  These are included to allow sophisticated users greater control over the optimization; the defaults seem to work well for most, but not all sequences. For \code{meth="L-BFGS-B"}, some parameters are constrained to be non-negative, which is useful paramters which cannot truly be negative, such as \code{vstep} (random walk) and \code{omega} (stasis model). \cr

  Initial estimates to start the optimization come from analytical solutions based on assuming equal sampling error across samples and evenly spaced samples in time (functions \code{mle.GRW}, \code{mle.URW} and \code{mle.Stasis}). 
  
}

\value{
  A \code{paleoTSfit} object.
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. }
\author{ Gene Hunt }

\note{
  Standard errors computed from the Hessian matrix are reasonably accurate for \code{mstep} and \code{theta}, but not as useful for the vstep and omega because of the asymmetry of the log-likelihood surfaces. 
}

\seealso{\code{\link{logL.GRW}}, \code{\link{fit3models}}, \code{\link{opt.RW.Mult}}, \code{\link{sim.GRW}} , \code{\link{as.paleoTSfit}}}

\examples{
 ## generate data for a directional sequence
 y <- sim.GRW(ns=30, ms=1, vs=1)
 plot(y)
 m.rw<- opt.GRW(y)
 m.rwu<- opt.URW(y)
 m.sta<- opt.Stasis(y)

 ## print log-likelihoods; easier to use function fit3models()
 compareModels(m.rw, m.rwu, m.sta)
}

\keyword{models}
\keyword{ts}