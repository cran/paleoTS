\name{sim.punc}
\alias{sim.punc}
\alias{sim.sgs}
\alias{sim.GRW.shift}

\title{ Simulate evolutionary time-series with changing dynamics }
\description{
  These functions simulate evolutionary time-series in which the dynamics are not constant over time. 
  These correspond to unsampled punctuations (\code{sim.punc}), sampled punctuations (\code{sim.sgs}), 
  and random walks with parameter values change over time.
}

\usage{
sim.punc(ns = c(10, 10), theta = c(0, 1), omega=rep(0,length(theta)), 
         nn = rep(30, sum(ns)), tt = 0:(sum(ns)-1), vp = 1)
sim.sgs(ns = c(20, 20, 20), theta = 0, omega = 1, ms = 1, vs = 0.1, 
        nn=rep(30, sum(ns)), tt = 0:(sum(ns)-1), vp = 1)
sim.GRW.shift(ns = c(10, 10), ms = c(0, 1), vs = c(0.5, 0.5), 
              nn = rep(30, sum(ns)), tt = 0:(sum(ns)-1), vp = 1)
}

\arguments{
  \item{ns}{ vector with the number of samples in each segment }
  \item{theta}{ evolutionary optimum, Stasis model }
  \item{omega}{ evolutionary variance, Stasis model }
  \item{nn}{ a vector of sample sizes per sample }
  \item{tt}{ vector of sample ages, increasing from oldest to youngest }
  \item{vp}{ within-sample phenotypic variance }
  \item{ms}{ mean step in GRW model for puncuated interval}
  \item{vs}{ step variance in GRW model for puncuated interval}
}
\details{
These three models are elaborations of the standard random walk, directional evolution and stasis models
to allow for heterogeneous evolutionary dynamics within a sequence.  These heterogeneous dynamics include 
sampled punctuations (\code{sim.sgs}), unsampled punctuations {\code{sim.punc}} and shifting random walk
dynamics (\code{sim.GRW.shift}).
}

\value{
  A \code{paleoTS} object.
}

\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology} \bold{ 34}:360-377. \cr
	 }

\author{ Gene Hunt }


\seealso{ \code{\link{opt.punc}}, \code{\link{opt.sgs}}, \code{\link{opt.GRW.shift}} } 
\examples{
# illustrate differences between sampled and unsampled punctuations
x1<- sim.sgs(ms=0.4, omega=0.5)
x2<- sim.punc(ns=c(30,30), theta=c(0, mean(x1$mm[41:60])), omega=c(0.5,0.5))
layout(1:2)
plot(x1, nse=2, col="red", main="Sampled Punctuation")
rect(x1$tt[20],min(x1$mm), x1$tt[40], max(x1$mm))
text(x1$tt[30], max(x1$mm), "sampled\ntransition", font=3, pos=1, cex=0.8)
plot(x2, nse=2, col="black", main="Unsampled Punctuation")
rect(x2$tt[30],min(x2$mm),x2$tt[31],max(x2$mm))
text(x2$tt[31], mean(x2$mm), "unsampled\ntransition", font=3, pos=4, cex=0.8)
}

\keyword{models}
\keyword{ts}