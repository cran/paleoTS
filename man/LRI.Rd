\name{LRI}
\alias{LRI}

\title{
Log-Rate, Log-Interval (LRI) method of Gingerich
}

\description{
Gingerich (1993) introduced a method that plots on log-log scale, the rate and interval for each pair of samples in an evolutionary sequence.  
On this plot, the slope is interpreted as an indicator of evolutionary mode (-1 for stasis, 0.5 for random walk, 0 for directional), and the intercept is interpreted as a measur of the rate of evolution over one generation.
}

\usage{
LRI(x, gen.per.t = 1e+06, draw = TRUE)
}

\arguments{
  \item{x}{a \code{paleoTS} object}
  \item{gen.per.t}{the number of generations per unit time (e.g., 1e6 for yearly generations and time in \code{x} is in Myr)}
  \item{draw}{logical, if \code{TRUE}, a plot is produced}
}

\details{
Following Gingerich (1993), a robust line is fit through the points by minimizing the sum of absolute deviations.
}

\value{
A named vector of three elements: \code{Intercept}, \code{slope} and \code{GenerationalRate}
}

\references{
Gingerich, P.D. 1993. Quantification and comparison of evolutionary rates. \emph{American Journal of Science} \bold{ 293-A}:453--478. 
Gingerich, P.D. 2009. Rates of evolution. \emph{Annual Review of Ecology Evolution and Systematics} \bold{ 40}:657--675. 
}

\author{
Gene Hunt
}

\note{
This method was important in attempts to disentangle evolutionary tempo and mode.  I view likelihood-based methods as more informative,
 and in particular the estimation of 'Generational Rates' using LRI is compromised by sampling error (see the example).
}



\seealso{
 \code{\link{fit3models}}
}

\examples{
	xFast<- sim.GRW(ns=20, ms=0.1, vs=0.2)   # fast evolution
	xSlow<- sim.Stasis(ns=20, theta=10, omega=0)   # strict stasis! rates are actually zero
	wFast<- LRI(xFast, draw=FALSE)
	wSlow<- LRI(xSlow, draw=FALSE)
	## LRI usually assigns faster generational rate to the Strict Stasis!
	print(wFast[3],4)
	print(wSlow[3],4)
}

\keyword{models}
\keyword{ts}