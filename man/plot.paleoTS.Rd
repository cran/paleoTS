\name{plot.paleoTS}
\alias{plot.paleoTS}
\title{ Plots paleoTS objects }
\description{
  Plots paleontological time series, showing trait means (with error bars) over time. 
}
\usage{
plot.paleoTS(x, nse = 1, yl = NULL, pch = NULL, pool = TRUE, true.means = FALSE, add = FALSE, ...)
}
\arguments{
  \item{x}{ a \code{paleoTS} object }
  \item{nse}{ width of error bars, in standard errors }
  \item{yl}{ range of trait values to plot }
  \item{pch}{ plotting character }
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
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{as.paleoTS}} }

\examples{
y <- sim.GRW(50, 0, 1)
plot(y)
ys <- sub.paleoTS(y, k=0.2)		# same sequence, sub-sampled
plot (ys, add=TRUE, col="red", pch=4)
}
\keyword{models}
\keyword{ts}
\keyword{hplot}