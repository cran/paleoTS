\name{as.paleoTSfit}
\alias{as.paleoTSfit}

\title{
Class for fit to paleontological time-series models
}

\description{
Function combines useful information summarizing model fit.
}

\usage{
as.paleoTSfit(logL, parameters, modelName, method, K, n, se)
}

\arguments{
  \item{logL}{
log-likelihood of model 
}
  \item{parameters}{
named vector with maximum-likelihood parameter estimates
}
  \item{modelName}{
name of the model
}
  \item{method}{
parameterization used: ancestor-descendant (AD) or Joint
}
  \item{K}{
number of parameters in the model
}
  \item{n}{
sample size
}
  \item{se}{
standard errors of parameters from the hessian (optional)
}
}

\details{
This function is used by all the model-fitting routines (opt.XXX, fitXXX) to create standardized output
}
\value{
A list with the above items}

\author{
Gene Hunt
}

\note{
This function is not likely to be called directly by the user.
}


\seealso{
\code{\link{as.paleoTS}}
}

\keyword{models}
\keyword{ts}