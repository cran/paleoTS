\name{fitModeShift}
\alias{fitModeShift}
\alias{opt.joint.RW.Stasis}
\alias{opt.AD.RW.Stasis}
\alias{opt.joint.Stasis.RW}
\alias{opt.AD.Stasis.RW}
\alias{logL.joint.URW.Stasis}
\alias{logL.joint.GRW.Stasis}
\alias{logL.joint.Stasis.URW}
\alias{logL.joint.Stasis.GRW}


\title{
Fit models in which start in Stasis, and then shift to a random walk (or vice versa)
}
\description{
This function fits to a paleoTS object a model that starts in Stasis and shifts to either a unbiased or general random walk (models URW or GRW), or vice versa.
}
\usage{
fitModeShift(y, minb = 7, pool = TRUE, order = c("Stasis-RW", "RW-Stasis"),
             rw.model = c("URW", "GRW"), method = c("Joint", "AD"), 
             silent = FALSE, hess = FALSE, ...)

opt.joint.RW.Stasis(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), 
                    pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.AD.RW.Stasis(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), 
                 pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.joint.Stasis.RW(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), 
                    pool=TRUE, meth="L-BFGS-B", hess=FALSE)
opt.AD.Stasis.RW(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), 
                 pool=TRUE, meth="L-BFGS-B", hess=FALSE)
logL.joint.URW.Stasis(p, y, gg)
logL.joint.GRW.Stasis(p, y, gg)
logL.joint.Stasis.URW(p, y, gg)
logL.joint.Stasis.GRW(p, y, gg)
}

\arguments{
  \item{y}{a \code{ paleoTS} object}
  \item{minb}{ the minimum number of samples to require before and after the mode shift}
  \item{pool}{ logical, if TRUE, sample variances are replaced by their pooled estimate across all samples}
  \item{order}{ either "Stasis-RW" or "RW-Stasis", setting the order of the two modes}
  \item{rw.model}{ either "GRW" or "URW", setting the random walk model to be unbaised (URW) or general/biased/directional (GRW)}
  \item{method}{ parameterization to use: see Details}
  \item{silent}{ logical, if TRUE, results are not printed}
  \item{hess}{ logical, if TRUE, standard errors are computed from the Hessian matrix}
  \item{\dots}{ further arguments, to be passed to optimization functions}
  \item{gg}{ vector of group membership, from \code{shift2gg}}
  \item{cl}{ control list for optimization}
  \item{meth}{ optimization method, can be "BFGS" or "L-BFGS-B"}
  \item{p}{ parameter vectors for log-likelihood functions}  
}

\details{
Only the function \code{fitStasisRW} is likely to be used by the user.  The other functions are called directly or indirectly by this one.

The \code{method} argument refers to different parameterizations of the model. See the documentation under \code{\link{opt.joint.GRW}} for more information about the differences between these parameterizations.  
}
\value{
  A \code{paleoTSfit} object, with additional elements:
  \item{all.logl}{ log-likelihoods for all tested partitions of the series into segments}
  \item{GG}{ matrix of indices of initial samples of each tested segment configuration; each column of \code{GG} corresponds to the elements of \code{all.logl}}
}
\references{
	Hunt, G., M. J. Hopkins, and S. L. Lidgard 2015. Simple versus complex models of trait evolution and stasis as a response to environmental change.  \emph{PNAS} \bold{ 112}:4885--4890. \cr
}
\author{
Gene Hunt
}

\seealso{
\code{\link{as.paleoTSfit}}, \code{\link{bootSimpleComplex}}, \code{\link{fit9models}}
}
\examples{
	x<- sim.GRW(ns=30)
	m.urw<- fitSimple(x, model="URW")
	m.cplx<- fitModeShift(x, order="Stasis-RW", rw.model="URW", method="Joint")
	compareModels(m.urw, m.cplx)
}

\keyword{ models }
\keyword{ ts }
