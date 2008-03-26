\name{std.paleoTS}
\alias{std.paleoTS}

\title{ Standardize paleontological time series data }
\description{
 This function standardizes a \code{paleoTS} object, converting the data to within-populations standard deviation units. 
}
\usage{
std.paleoTS(y, zero = "start")
}
\arguments{
  \item{y}{ a paleontological time series }
  \item{zero}{ if \code{zero = "start"}, trait values are translated such that the initial starting value of the sequence is set to be zero.  All other values are ignored. }
}
\details{
  The standardization expresses each sample mean as the deviation from the overall mean, divided by the pooled within-sample variance.  Sample variances are also divided by the pooled sample variance. \cr
  Essentially, this converts paleontological time-series data into standard deviation units, similar to the computation of evolutionary rates in haldanes.  This operation \emph{does not} change the relative fit of models, but it does facilitate the comparison of parameter estimates across time-series with different units. 
}
\value{
  the converted \code{paleoTS} object.
}
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. 
}

\author{ Gene Hunt } 

\seealso{ \code{\link{ln.paleoTS}} }

\examples{
 y<- sim.GRW(20, 0, 1, vp=10)
 ys<- std.paleoTS(y)
 fit3models(y)
 fit3models(ys)  # note Akaike weights do not change from standardization
}

\keyword{models}
\keyword{ts}