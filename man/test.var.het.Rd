\name{test.var.het}
\alias{test.var.het}
\alias{pool.var}

\title{ Variance heterogeneity test }
\description{
  Tests for heterogeneity of sample variances in a \code{paleoTS} object.
}
\usage{
test.var.het(y, method = "Bartlett")
pool.var(y, nn=NULL, ret.paleoTS=FALSE)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{method}{ test to be used; currently only Bartlett's test is implemented   }
  \item{nn}{ if y is a vector of variances, nn is a vector of sample sizes }
  \item{ret.paleoTS}{ logical, indicating if the function should return a new \code{paleoTS} object with variances replaced by pooled variance.}
}
\details{
  Tests for variance heterogeneity among samples.  In the absence of evidence for heterogeniety, it may be desirable to replace individual estimates of phenotypic variance with a variance estimate pooled over all samples. \cr
  Function \code{pool.var} is used internally in \code{test.var.het}.  It returns the value of the pooled variance. 

}
\value{
Function \code{pool.var} either returns a \code{paeloTS} object, or the pooled variance. \cr  
Function \code{test.var.het} returns a list with the following elements relevant to Bartlett's test
  \item{stat }{ test statistic for Bartlett's test}
  \item{p.value }{ P-vlaue for statistical test}
  \item{df}{ degrees of freedom, equal to one fewer than the number of samples}
}
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601. 
\cr
Sokal, R. and F. J. Rohlf 1995.  Biometry, 3rd Ed.  }

\author{ Gene Hunt } 


\examples{
# generate simulated data with no real variance heterogeneity
y <- sim.GRW(20, 0, 1)
bart <- test.var.het(y)   # nonsignificant 95% of the time!
ys<- pool.var(y, ret.paleoTS=TRUE)
cat (ys$vv)		# note all variances now the same
}
\keyword{models}
\keyword{ts}