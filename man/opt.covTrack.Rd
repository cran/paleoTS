\name{opt.covTrack}
\alias{opt.covTrack}
\alias{opt.covTrack.Mult}
\alias{logL.covTrack}
\alias{logL.Mult.covTrack}

\title{
Covariate-tracking model
}

\description{
Functions to fit a model in which a phenotpyic trait tracks changes in a measured covariate (e.g., body size tracks temperature changes).
}
\usage{
opt.covTrack(y, z, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.covTrack.Mult(yl, zl, cl = list(fnscale = -1), pool = TRUE, hess = FALSE)
logL.covTrack(p, y, z)
logL.Mult.covTrack(p, yl, zl)
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
These functions model the changes in a phenotpyic trait as linear functions of changes in a covariate, \code{z}.
Each change in z (=dz) yields an expected change in the trait equal to \emph{b*dz}, with normally distributed residual variance \emph{evar}.
If \code{z} is equal in length to the trait sequence, first differences are taken.  If \code{z} is of length one less than the trait sequence, it is assumed that the user has already taken the appropriate difference.

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

\note{
These functions use the "AD" parameterization.
}



\seealso{
\code{\link{as.paleoTSfit}}, \code{\link{opt.GRW}}
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