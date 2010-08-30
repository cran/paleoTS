\name{sim.covTrack}
\alias{sim.covTrack}

\title{ Simulate time-series that tracks a covariate 
}

\description{
Function to simulate a time-series of trait values in which each evolutionary increment is linearly related to changes in a covariate
}

\usage{
sim.covTrack(ns = 20, b = 1, evar = 0.1, z, nn = rep(20, times = ns), tt = 0:(ns-1), vp = 1)
}

\arguments{
  \item{ns}{ number of samples in time-series }
  \item{b}{ slope of the relationship between evolutionary changes and changes in the covaiate}
  \item{evar}{ residual variance around covariate-trait relationship}
  \item{z}{ a measured covariate }
  \item{nn}{ vector of the number of individuals in each sample }
  \item{tt}{ vector of sample ages, increases from oldest to youngest }
 \item{vp}{ within-population trait variance }
}

\details{
See \code{opt.covTrack} for a description of the model.
}

\value{
An object of class \code{paleoTSfit}}

\references{
Hunt, et al. 2010. Climate-driven body-size trends in the ostracod fauna of the deep Indian Ocean, \emph{Palaeontology} \bold{ 53}:1255--1268. 
}

\author{
Gene Hunt
}


\seealso{
\code{\link{opt.covTrack}}, \code{\link{as.paleoTSfit}}
}
\examples{
 z<- rnorm(20)
 x<- sim.covTrack(ns=20, b=2, evar=0.3, z)
 plot(diff(z), diff(x$mm), xlab="Change in covariate", ylab="Change in Trait")
 abline(h=0, lty=3)
 abline(v=0, lty=3)

}

\keyword{models}
\keyword{ts}