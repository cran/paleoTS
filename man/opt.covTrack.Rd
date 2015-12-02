\name{opt.covTrack}
\alias{opt.covTrack}
\alias{opt.covTrack.Mult}
\alias{logL.covTrack}
\alias{logL.Mult.covTrack} 
\alias{opt.joint.covTrack} 
\alias{opt.joint.covTrack.Mult} 
\alias{logL.joint.covTrack} 
\alias{logL.Mult.joint.covTrack} 

\title{
Covariate-tracking model
}

\description{
Functions to fit a model in which a phenotpyic trait tracks changes in a measured covariate (e.g., body size tracks temperature changes).
}
\usage{
opt.covTrack(y, z, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.joint.covTrack(y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
opt.covTrack.Mult(yl, zl, cl = list(fnscale = -1), pool = TRUE, hess = FALSE)
opt.joint.covTrack.Mult(yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)

logL.covTrack(p, y, z)
logL.joint.covTrack(p, y, z)
logL.Mult.covTrack(p, yl, zl)
logL.Mult.joint.covTrack(p, yl, zl)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{z}{ a measured covariate.  See \code{Details} about its required length.}
  \item{pool}{ logical, if \code{TRUE}, variance are pooled across samples }
  \item{cl}{ control list, passed to function \code{optim} }
  \item{meth}{ optimization method, passed to function \code{optim} }
  \item{hess}{ logical, indicating whether to calculate standard errors from the Hessian matrix }
  \item{yl}{ a \code{list} of paleoTS objects } 
  \item{zl}{ a \code{list} of covariates, corresponding to \code{yl} }
  \item{p}{ a vector of parameter values } 
}

\details{
These functions fit a model in trait evolution tracks a covariate (\emph{z}) over time.  There are two parameterizations: "Joint" and "AD"; all functions without "joint" in their names use the "AD" parameterization.  The "joint" parameterization assumes the trait values are a linear function of the covariate, whereas the "AD" parameterization assumes that \emph{changes} in the traits are a linear function of \emph{changes} in the covariate.

For the "AD" parameterization, \code{z} can be equal in length to y, or have one fewer element.  If the former, first differences are taken.  If \code{z} is of length one less than the trait sequence, it is assumed that the user has already taken the appropriate difference.  For the "Joint" parameterization, \code{z} must be of the same length as \code{y}. 

For the \code{Mult} versions of these, the function estimates the model assuming the \emph{same} relationship exists across all sequences.
}

\value{
An object of class \code{paleoTSfit}
}

\references{
Hunt, et al. 2010. Climate-driven body-size trends in the ostracod fauna of the deep Indian Ocean, \emph{Palaeontology} \bold{ 53}:1255--1268. 
}

\author{
Gene Hunt
}

\section{Note}{
	The convenience function \code{fitSimple} can also be called to fit this model.	
}

\section{Warning}{
The "Joint" parameterization does not de-trend the series in any way, and so is vulnerable to type I error if both \code{y} and \code{z} are trended.  If so, an warning message is given.
}


\seealso{
\code{\link{as.paleoTSfit}}, \code{\link{opt.GRW}}, \code{\link{fitSimple}}
}

\examples{
 z<- rnorm(20)
 x<- sim.covTrack(ns=20, b=2, evar=0.1, z=z)
 plot(diff(z), diff(x$mm), xlab="Change in covariate", ylab="Change in Trait")
 abline(h=0, lty=3)
 abline(v=0, lty=3)
 mct<- opt.covTrack(x, z)
 print(round(mct$parameters,2))  # should be close to generating values
 print(x$genpar)
}

\keyword{models}
\keyword{ts}