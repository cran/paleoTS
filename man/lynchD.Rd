\name{lynchD}
\alias{lynchD}

\title{ Compute rate metric from Lynch (1990) }
\description{
 This function computes D, the rate metric proposed by Lynch (1990).  This metric derives from the random walk model, with \code{D = Vstep/Vp}, where \code{Vstep} is the
 step variance of the unbiased random walk, and \code{Vp} is the within sample variance, pooled among samples.  Under mutation - drift equilibrium, \code{D} is expected to range
 approximately between 5e-5 and 5e-3. 
}
\usage{
lynchD(y, gen.per.t = 1e+06, pool = TRUE, ...)
}
\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{gen.per.t}{ the number of generations per unit time }
  \item{pool}{ logical, if the variance should be pooled across samples in estimating \code{Vstep} }
  \item{\dots}{ further arguments, passed to \code{opt.URW}}
}
\details{
 The \code{gen.per.t} argument indicates the number of organismal generations for each unit of time with which the time-series \code{y} was measured.  For example, is \code{y$tt} is measured in millions of years, and the species has annual generations, \code{gen.per.t} would be 1e6.
\cr
 Estimation of \code{Vstep} is done on the original time scale of \code{y}, and then the values are converted to generational time scales using \code{gen.per.t}.  This is for numerical reasons, as it avoids computations on possible very low numbers. 
}
\value{
  \item{D}{  value of rate metric }
  \item{pooled.var}{ value of pooled within-sample variance }
  \item{gen.per.t}{ number of generations per unit time }
  \item{vstep} { computed \code{Vstep}, at the original time scale of \code{y}}
  \item{drift.range}{ exepcted minimum and maximum values of \code{D} consistent with neutral evoltuion }
  \item{result}{ conclusion reached about the plausibility of neutral evolution }
}

\references{ Lynch, M. 1990. The rate of morphological evolution in mammals from the standpoint of the neutral expectation. \emph{The American Naturalist} \bold{136}:727--741. }

\author{ Gene Hunt }

\seealso{ \code{\link{opt.URW}} }

\examples{
 y<- sim.GRW(ns=20, ms=0, vs=1e-4, tt=seq(0,1e6, length.out=20))  # per-year simulation
 lynchD(y, gen.per.t=1)		# 1 generation per year
}	
\keyword{models}
\keyword{ts}