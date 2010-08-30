\name{plot.paleoTS}
\alias{plot.paleoTS}
\title{ Plots paleoTS objects }
\description{
  Plots paleontological time series, showing trait means (with error bars) over time. 
}
\usage{
\method{plot}{paleoTS}(x, nse = 1, pt.ch = 19, yl = NULL, pool = FALSE, true.means = FALSE, add = FALSE, ...)
}
\arguments{
  \item{x}{ a \code{paleoTS} object }
  \item{nse}{ size of error bars, in standard errors }
  \item{pt.ch}{ plotting character }
  \item{yl}{ range of trait values to plot }
  \item{pool}{ logical, indicating if should pool variances among samples in computing standard errors }
  \item{true.means}{ logical indicating if true means (from simulated data) should also be plotted }
  \item{add}{ logical, if \code{TRUE} data are added to existing plot }
  \item{\dots}{ other arguments sent to \code{plot} }
}
\details{
 True means are only available when the data are simulated.  If plotted, they are indicated by a dashed blue line.  
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
plot (ys, add=TRUE, col="red", pt.ch=4)
}
\keyword{models}
\keyword{ts}
\keyword{hplot}