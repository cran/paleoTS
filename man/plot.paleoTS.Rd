\name{plot.paleoTS}
\alias{plot.paleoTS}

\title{ Plots paleoTS objects }
\description{
  Plots paleontological time series, showing trait means (with error bars) over time. 
}
\usage{
\method{plot}{paleoTS}(x, nse = 1, pool = FALSE, add = FALSE, modelFit=NULL, pch=19, lwd=1.5, ylim=NULL, ...)
}
\arguments{
  \item{x}{ a \code{paleoTS} object }
  \item{nse}{ size of error bars, in standard errors }
  \item{pool}{ logical, indicating if should pool variances among samples in computing standard errors }
  \item{add}{ logical, if \code{TRUE} data are added to existing plot }
  \item{modelFit}{ optional \code{paleoTSfit} object with model fit results to add to plot }
  \item{pch}{plotting symbol}
  \item{lwd}{line width for plotting}
  \item{ylim}{y limit for plotting (optional)}
  \item{\dots}{ other arguments sent to \code{plot} }
}
\details{
 If \code{modelFit} is not NULL, then the expectations of that model, plus 95% probability intervals, are added to the plot.    
}

\value{
 No values are returned, the results are plotted. 
}

\author{ Gene Hunt }

\seealso{ \code{\link{as.paleoTS}} }

\examples{
y <- sim.GRW(50, 0, 1)
plot(y)
ys <- sub.paleoTS(y, k=0.2)		# same sequence, sub-sampled
plot (ys, add=TRUE, col="red")
}
\keyword{models}
\keyword{ts}
\keyword{hplot}