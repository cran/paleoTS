% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basics.R
\name{std.paleoTS}
\alias{std.paleoTS}
\title{Convert time-series to standard deviation units}
\usage{
std.paleoTS(y, center = c("mean", "start"))
}
\arguments{
\item{y}{a \code{paleoTS} object}

\item{center}{optional translation of time-series according to "mean" or
"start"; see Details}
}
\value{
the converted \code{paleoTS} object
}
\description{
Convert time-series to standard deviation units
}
\details{
The standardization expresses each sample mean as the deviation
  from the overall mean, divided by the pooled within-sample standard deviation. Sample
  variances are also divided by the pooled sample variance. \cr \cr Essentially,
  this converts paleontological time-series data into standard deviation
  units, similar to the computation of evolutionary rates in haldanes.  This
  operation \emph{does not} change the relative fit of models, but it does
  facilitate the comparison of parameter estimates across time-series of
  traits measured in different units. \cr \cr If argument \code{center} = "start"
  the time-series is translated such that the trait mean of the first sample
  is zero.
}
\examples{
x <- sim.Stasis(ns = 8, theta = 1, omega = 4, vp = 2)
xs <- std.paleoTS(x, center = "start")
plot(x, ylim = range(c(x$mm, xs$mm)))
plot(xs, col = "red", add = TRUE)
legend(x = "topright", c("unstandardized", "standardized"), lty=1, col=c("black", "red"), cex=0.7)
}
