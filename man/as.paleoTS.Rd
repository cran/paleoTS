\name{as.paleoTS}
\alias{as.paleoTS}
\alias{read.paleoTS}

\title{ Paleontological time-series class }
\description{
  A class for storing information on evolutionary time-series data.
}
\usage{
as.paleoTS(mm, vv, nn, tt, MM = NULL, genpars = NULL, label = NULL, start.age = NULL, oldest=c("first", "last"), reset.time=TRUE)
read.paleoTS(file = NULL, oldest = "first", reset.time=TRUE, ...)
}

\arguments{
  \item{mm}{ vector of sample means }
  \item{vv}{ vector of sample variances }
  \item{nn}{ vector of sample sizes }
  \item{tt}{ vector of sample ages }
  \item{MM}{ vector of true population means (for simulated time series) }
  \item{genpars}{ generating parameters (for simulated time series) }
  \item{label}{ optional label describing time series data }
  \item{start.age}{ optional, indicating age of oldest sample (e.g., in millions of years ago) }
  \item{file}{ file with \code{paleoTS} data; if \code{NULL}, then opens \code{file.choose} dialog }
  \item{oldest}{ by default, samples are listed from oldest to youngest.  If file is in reverse order, use \code{oldest = "last"}}
  \item{reset.time}{logical; if TRUE, then change time scale to start at t=0 and adjust \code{start.age} accordingly.}
  \item{\dots}{ further arguments passed to \code{read.table}}
}
\details{
   Function \code{as.paleoTS()} combines time series data into an object of class \code{paleoTS}.  This function will usually not be used directly; \code{read.paleoTS()} is more convenient for getting the relevant data from text files.
\cr
  If no file name is given for \code{read.paleoTS}, the user will be prompted to select a file using the \code{file.choose()} interactive prompt.  Samples are assumed to be ordered from oldest to youngest (if not, use \code{oldest="last"}), with ages indicating time elapsed from the beginning of the sequence, rather than geological age before the present.  If sample ages decrease through the sequence, as if given in millions of years ago, \code{tt} will automatically be converted to time elapsed from the beginning of the sequence. Generally, setting \code{reset.time} = TRUE is good practice.
}

\note{
  The analysis estimates the contribution of sampling noise to the observed differences between samples.  It does this assuming that one is analyzing sample means, which have sampling variances equal to variance divided by sample size, \code{vv/nn}.  If one is interested in analyzing statistics other than the sample mean (medians, quantiles, or other statistics), use the the following procedure: set the statistic in question as the \code{mm} values, replace \code{vv} with a vector of the squared standard errors for each estimate, and set all values of \code{nn} to one.   
}

\value{
  \item{mm}{ vector of sample means }
  \item{vv}{ vector of sample variances }
  \item{nn}{ vector of sample sizes }
  \item{tt}{ vector of sample ages }
  \item{MM}{ vector of true population means (for simulated time series) }
  \item{genpars}{ generating parameters (for simulated time series) }
  \item{label}{ optional label describing time series data }
  \item{start.age}{ optional, indicating age of oldest sample (e.g., in millions of years ago) }
}

\section{Warning}{
Some functions do assume that \code{reset.time=TRUE}, that is that time (\code{tt}) starts at zero and increases from there.  Do not change this argument unless you know what you are doing. 		
}

\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{plot.paleoTS}}  }
\examples{
## create paleoTS object
 y <- as.paleoTS(mm=rnorm(10), vv=rep(1,10), nn=rep(25,10), tt=1:10, label="white noise time series")
 plot(y)
}

\keyword{models}
\keyword{ts}