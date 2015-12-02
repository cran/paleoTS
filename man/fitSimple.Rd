\name{fitSimple}
\alias{fitSimple}


\title{
Fit simple models of trait evolution
}
\description{
A convenience wrapper function that can calls other \code{paleoTS} functions to fit the unbiased random walk (URW), general random walk (GRW), Stasis, Strict Stasis, Ornstein-Uhlenbeck (OU) and covariate-tracking (covTrack) models.
}
\usage{
fitSimple(y, model = c("GRW", "URW", "Stasis", "StrictStasis", "OU", "covTrack"),
          method = c("Joint", "AD"), pool = TRUE, z = NULL, hess = FALSE)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{model}{ the model to be fit; one of c("GRW", "URW", "Stasis", "OU", "covTrack")}
  \item{method}{ parameterization to use: \code{Joint} or \code{AD}; see Details}
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{z}{ the covariate variable; only used for the covTrack model}
  \item{hess}{ logical, indicating whether to calculate standard errors from the Hessian matrix }
   
}

\details{
For the covariate-tracking model, z should be a vector of length \emph{n} when \code{method="Joint"} and \emph{n}-1 when \code{method="AD"}, where \emph{n} is the number of samples in \code{y}.

Note that the AD method has not been implemented for the OU model.  The Joint method seems to do rather better for this model, anyway.

\code{fitMult} fits these models (excpet for the OU model) over multiple \code{paleoTS} objects, under the assumption that the same model applies to all the trait sequences.  Parameters other than the stasis mean (\code{theta}) and the ancestral state (\code{anc}) are also assumed to be shared among sequences.
}

\value{
  A \code{paleoTSfit} object.
}

\references{
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. \cr
	Hunt, G., M. J. Hopkins, and S. L. Lidgard 2015. Simple versus complex models of trait evolution and stasis as a response to environmental change.  \emph{PNAS} \bold{ 112}:4885--4890. \cr
}

\author{
Gene Hunt
}

\seealso{
\code{\link{opt.GRW}}, \code{\link{opt.joint.GRW}}, \code{\link{opt.covTrack}}, \code{\link{opt.joint.GRW}} 
}

\examples{
	x<- sim.Stasis(ns=30, theta=10, omega=1)
	s1<- fitSimple(x, model="URW")
	s2<- fitSimple(x, model="Stasis")
	s3<- fitSimple(x, model="StrictStasis")
	compareModels(s1, s2, s3)
}


\keyword{ models }
\keyword{ ts }
