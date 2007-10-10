\name{ln.paleoTS}
\alias{ln.paleoTS}

\title{ Log transform paleontological time series data }
\description{
 This function performs an approximate log-transformation (base e) of a \code{paleoTS} object.
}
\usage{
ln.paleoTS(y)
}
\arguments{
  \item{y}{ a paleontological time series }
}
\details{
  For a random variable x, its approximate mean on a log scale is the log of its untransformed mean.  The approximate variance on a log scale is equal to the squared coefficient of variation.
}
\value{
  the converted \code{paleoTS} object.
}
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. 
\cr
Lewontin, R. 1966. On the measurement of relative variability.  \emph{Systematic Zoology} \bold{15}:141--142. }
}

\author{ Gene Hunt } 

\seealso{ \code{\link{std.paleoTS}} }

\examples{
 y <- sim.GRW(20, 20, 1)
 yl <- ln.paleoTS(y)
 print (yl)
}

\keyword{models}
\keyword{ts}