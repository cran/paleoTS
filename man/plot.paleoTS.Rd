\name{plot.paleoTS}
\alias{plot.paleoTS}
\title{ Plots paleoTS objects }
\description{
  Plots paleontological time series, showing trait means (with error bars) over time. 
}
\usage{
plot.paleoTS(x, nse = 1, yl = NULL, pch = 19, true.means = FALSE, start.time = NULL, add = FALSE, ...)
}
\arguments{
  \item{x}{ a \code{paleoTS} object }
  \item{nse}{ width of error bars, in standard errors }
  \item{yl}{ range of trait values to plot }
  \item{pch}{ plotting character }
  \item{true.means}{ logical indicating if true means (from simulated data) should also be plotted }
  \item{start.time}{ Value indicating the absolute age of the oldest sample  }
  \item{add}{ logical, if \code{TRUE} data are added to existing plot }
  \item{\dots}{ other arguments sent to \code{plot} }
}
\details{
 True means are only available when the data are simulated.  If plotted, they are indicated by a dashed blue line.  
}
\value{
 No values are returned, the results are plotted. 
}
\references{ Hunt, G. 2007. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{as.paleoTS}} }

\examples{
y <- sim.rw(20, 0, 1)
plot(y, nse=2, true.means=TRUE, start.time=100)
}
\keyword{models}
\keyword{ts}
\keyword{hplot}