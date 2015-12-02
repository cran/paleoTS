\name{modelCurves}
\alias{modelCurves}

\title{
Function computes model expectations and 95% probability interval for some models.
}
\description{
This function is used by \code{plot.paleoTS} to add the model expectations and probability intervals for (some) models.  It is not likely to be called directly by the user.
}

\usage{
modelCurves(x, w, np = 500)
}

\arguments{
  \item{x}{a \code{paleoTS} object}
  \item{w}{fitted model results in the form of a \code{paleoTSfit} object}
  \item{np}{the number of points to use for plotting the expectation}
}

\details{
The function determines what to draw from the \code{modelName} of \code{w}, along with its parameter values.
}

\value{
  A list with the following vectors: \code{tt}, the times plotted; \code{ee}, the model expectations; \code{ll} and \code{uu}, the lower and upper 95 percent probability limits for the model.
}

\author{
Gene Hunt
}

\note{
If the \code{modelName} of \code{w} is not recognized, \code{ee} is set to NA and nothing is added to the plot.  At present, the following models are implemented: GRW, URW, Stasis, unsampled punctuations (as by \code{fitGpunc}).
}

\examples{
 x1<- sim.punc(theta=c(0,10), omega=c(2,2))
 w1<- fitGpunc(x1, ng=2)
 plot(x1, modelFit=w1) 
 
 x2<- sim.GRW(ns=20, ms=1, vs=0.4)
 w2<- opt.GRW(x2)
 plot(x2, modelFit=w2)

}

\keyword{models}
\keyword{ts}% 
