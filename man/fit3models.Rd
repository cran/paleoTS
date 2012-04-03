\name{fit3models}
\alias{fit3models}

\title{Do model fits for three standard evolutionary models}

\description{
  Functions fits three models to an evolutionary time series: (1) general random walk (=directional evolution), (2) unbiased random walk, and (3) stasis. 
}

\usage{
fit3models(y, silent = FALSE, method=c("AD", "Joint"), ...)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{silent}{ logical, if TRUE, results are not printed }
  \item{method}{ parameterization to use: see Details}
  \item{\dots}{ further arguments, passed to optimization functions}
}

\details{
The \code{method} argument refers to different parameterizations of the model. See the documentation under \code{\link{opt.joint.GRW}} for more information about the differences between these parameterizations.  
}

\note{
For the \dots argument, only the arguments \code{pool} and \code{hess} will work; specifying \code{meth} in particular will cause an error.  If finer control is needed, the optimization functions should be used individually.
}

\value{
  If \code{silent=FALSE}, function \code{fit3models} prints and returns a dataframe with log-likelihoods, numbers of parameters, AICc scores and Akaike weights for the three models.
  If \code{silent=TRUE}, a list is returned with element 'modelFits' that has the same dataframe, plus an element 'parameters' with the sub-elements for all the parameter estimates.
}

\references{ 
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. \cr
	Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based assessment of modes, rates, punctuations and process..  \emph{In} R.K. Bambach and P.H. Kelley, eds.  From Evolution to Geobiology: Research Questions Driving Paleontology at the Start of a New Century:578--601. \cr
	}


\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}}, \code{\link{opt.joint.GRW}} }

\examples{
 ## show difference in parameterizations
 ### example 1, sequence with a strong trend ###
 # two parameterizations usually yield similar Akaike weights under these conditions
 x1<- sim.GRW(ns=10, ms=1, vs=0.5)
 res1AD<- fit3models(x1, method='AD')
 res1Joint<- fit3models(x1, method='Joint')

 ## example 2, longer & noisy directional walk
 ## joint parameterization often is better at correctly favoring GRW under these conditions
 x2<- sim.GRW(ns=20, ms=0.2, vs=0.1)  # step variance relatively low compared to sampling error == Noisy
 res2AD<- fit3models(x2, method='AD')
 res2Joint<- fit3models(x2, method='Joint')


}
\keyword{models}
\keyword{ts}