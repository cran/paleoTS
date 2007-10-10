\name{sub.paleoTS}
\alias{sub.paleoTS}

\title{ Subset an evolutionary time series }
\description{
  This function subsets a \code{paleoTS} object, returning a shorter series of specified or random samples from the original.
}
\usage{
sub.paleoTS(y, ok = NULL, k = 0.1)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{ok}{ vector indicating which samples to retain; can be logical or numeric }
  \item{k}{ proportion of samples to retain, chosen randomly }
}
\details{
  If \code{ok = NULL}, the subsetting is done by retaining a proportion \code{k} of the initial samples, chosen at random.  If \code{ok} is specified, those samples indicated will be retained.
}
\value{
 the subetted \code{paleoTS} object
}
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601.  }

\author{ Gene Hunt }

\seealso{ \code{\link{as.paleoTS}} }

\examples{
 y <- sim.GRW(ns=100, ms=0, vs=0.3)
 plot(y)
 ys1<- sub.paleoTS(y, ok=y$mm > mean(y$mm))  # subsetting with logical ok
 ys2<- sub.paleoTS(y, ok=1:10)				 # subsetting with numeric ok
 ys3<- sub.paleoTS(y, k=0.1)				 # subsetting random 10%
 plot(ys1, add=TRUE, col="red")
 plot(ys2, add=TRUE, col="blue")
 plot(ys3, add=TRUE, col="grey", cex=2, pch=1)
}

\keyword{models}
\keyword{ts}