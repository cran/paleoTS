\name{as.paleoTS}
\alias{as.paleoTS}
\alias{read.paleoTS}

\title{ Paleontological time-series class }
\description{
  A class for storing information on evolutionary time-series data.
}
\usage{
as.paleoTS(mm, vv, nn, tt, MM = NULL, genpars = NULL, label = "", start.age = NULL)
read.paleoTS(file = NULL, hh = FALSE, oldest = "first", ...)
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
  \item{hh}{ logical; if TRUE, file with paleoTS data has a header}
  \item{oldest}{ by default, samples are listed from oldest to youngest.  If file is in reverse order, use \code{oldest = "last"}}
  \item{\dots}{ further arguments passed to \code{read.table}}
}
\details{
   Function \code{as.paleoTS()} combines time series data into an object of class \code{paleoTS}.  This function will usually not be used directly; \code{read.paleoTS()} is used to get the relevant data from text files, and \code{sim.GRW()} and \code{sim.Stasis()} are used to simulate paleontological time-series from generating parameters.
\cr
  If no file name is given for \code{read.paleoTS}, the user will be prompted to select a file using the \code{file.choose()} interactive prompt.  Samples are assumed to be ordered from oldest to youngest (if not, use \code{oldest="last"}), with ages indicating time elapsed from the beginning of the sequence, rather than geological age before the present.  If sample ages decrease through the sequence, as if given in millions of years ago, \code{tt} will automatically be converted to time elapsed from the beginning of the sequence. 
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
\references{ Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. }
\author{ Gene Hunt }

\seealso{ \code{\link{plot.paleoTS}}  }
\examples{
## create paleoTS object
 y <- as.paleoTS(mm=rnorm(10), vv=rep(1,10), nn=rep(25,10), tt=1:10, label="white noise time series")
 plot(y)
}

\keyword{models}
\keyword{ts}