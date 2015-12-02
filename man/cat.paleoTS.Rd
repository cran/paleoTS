\name{cat.paleoTS}
\alias{cat.paleoTS}
\alias{split4punc}
\alias{shifts}
\alias{shift2gg}

\title{ Miscellaneous functions used internally for punctuations }

\description{
These functions are used internally in the simulation and fitting of models with heterogeneous evolutionary 
dynamics.  They generally will not be used directly by users.
}

\usage{
cat.paleoTS(y)
split4punc(y, gg, overlap=TRUE)
shifts(ns, ng, minb=7)
shift2gg(ss, ns)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{gg}{ a vector of integers from 1 to the number of segments, indicating for each sample, the segment to which it belongs}
  \item{overlap}{ if TRUE, the initial sample of each segement starts with the last sample of the previous segment }
  \item{ns}{ the number of samples in the entire sequence}
  \item{ng}{ the number of 'groups', i.e., independently evolving evolutionary segments}
  \item{minb}{ the minimum number of samples to include in a segment}
  \item{ss}{ an integer, or vector of integers, indicating the index of the initial sample of segments }
}
\details{
  Function \code{cat.paleoTS} concatenates a list of multiple \code{paleoTS} objects into a single time-series.  
  Note that the ages (\code{tt}) of the resulting time-series are read directly from the constituent time-series. \cr
  
  Function \code{split4punc} divides a \code{paleoTS} time-series into separate component time-series according to the grouping vector \code{gg}. \cr
  
  Function \code{shifts} returns a matrix of all allowable shift points for dividing \code{ns} samples into \code{ng} segments,  subject to the constraint
  that each segement has at least \code{minb} samples.  Function \code{shifts2gg} converts a vector of shift points into a grouping vector.
}

\value{
 Function \code{cat.paleoTS} returns a single \code{paleoTS} object, and \code{split4punc} returns a list of multiple \code{paleoTS} objects. \cr
 
 Function \code{shifts} returns a matrix of shift points--the index of samples that initiate new segments.  This matrix has \code{ng}-1 rows, where the ith row 
 corresponds to the index of the initial sample for segment i+1.  Each column represents one of the possible shift point combinations. Function \code{shift2gg} 
 returns a grouping vector from these initial shift points.
}

\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology}\bold{32}:578--601.  \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology}\bold{ 34}:360-377.  \cr
}
\author{ Gene Hunt }

\seealso{ \code{\link{fitGpunc}}, \code{\link{sim.punc}}, \code{\link{opt.GRW.shift}} }

\examples{
 # show all possible shift points for 20 samples divided into two semgments 
 # (assuming each has at least 6 samples)
 GG<- shifts(ns=20, ng=2, minb=6)
 print (GG)
 # show resulting grouping vector from first of these
 gg<- shift2gg(GG[,1], ns=20)
 print (gg)
 
}
\keyword{models}
\keyword{ts}
