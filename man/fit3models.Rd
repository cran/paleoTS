\name{fit3models}
\alias{fit3models}
\alias{fit3models.joint}

\title{Do model fits for three standard evolutionary models}

\description{
  Functions fits three models to an evolutionary time series: (1) general random walk (=directional evolution), (2) unbiased random walk, and (3) stasis. 
}

\usage{
fit3models(y, pool = TRUE, silent = FALSE, wts = "AICc")
fit3models.joint(y, pool = TRUE, silent = FALSE, wts = "AICc")
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{pool}{ logical, if TRUE, variances are pooled across samples }
  \item{silent}{ logical, if TRUE, results are not printed }
  \item{wts}{ what version of AIC to use for calculation of Akaike weights; either \code{AIC} or \code{AICc} }
}

\details{
The two versions refer to different parameterizations of the model.  The plain version, \code{fit3models}, works from evolutionary differences between 
adjacent samples, and calls the functions \code{opt.GRW}, \code{opt.URW} and \code{opt.Stasis}.  The joint version, \code{fit3models.joint}, considers all samples jointly rather than 
working with the differences between samples, and calls the functions \code{opt.joint.GRW}, \code{opt.joint.URW} and \code{opt.joint.Stasis}. \cr
See the documentation under \code{\link{opt.joint.GRW}} for more information about the differences between these parameterizations.  
}
\value{
  If \code{silent=FALSE}, function \code{fit3models} returns a vector of parameter estimates; if \code{silent=TRUE}, a list with the following components is returned, with models listed in the order of general random walk, unbiased random walk, and stasis: 
  \item{aic }{Akaike information criterion}
  \item{aicc }{AIC modified for low sample sizes}
  \item{logl}{log-likelihoods of the three models}
  \item{hats}{parameter estimates for general random walk (mstep, vstep), unbiased random walk (vstep) and stasis (theta, omega) models}
  \item{ak.wts}{vector of Akaike weights}
}

\references{ 
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. \cr
	Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{62}:700--710. \cr   
	}


\author{ Gene Hunt }

\seealso{ \code{\link{opt.GRW}}, \code{\link{opt.joint.GRW}} }

\examples{
 ## show difference in parameterizations
 ### example 1, sequence with a strong trend ###
 # two parameterizations usually yield similar Akaike weights under these conditions
 x1<- sim.GRW(ns=10, ms=1, vs=0.5)
 plot(x1, nse=2)
 res1<- fit3models(x1, silent=TRUE)
 res1a<- fit3models.joint(x1, silent=TRUE)
 res1.tab<- rbind(res1$ak.wts, res1a$ak.wts)
 row.names(res1.tab)<- c("fit3models", "fit3models.joint")
 cat ("\n-- Akaike Weights from different parameterizations --\n")
 print (res1.tab)

 ## example 2, longer & noisy random walk
 ## joint parameterization often is better at correctly favoring URW under these conditions
 x2<- sim.GRW(ns=20, ms=0, vs=0.1)  # step variance relatively low compared to sampling error == Noisy
 plot(x2, nse=2)
 res2<- fit3models(x2, silent=TRUE)
 res2a<- fit3models.joint(x2, silent=TRUE)
 res2.tab<- rbind(res2$ak.wts, res2a$ak.wts)
 row.names(res2.tab)<- c("fit3models", "fit3models.joint")
 cat ("\n-- Akaike Weights from different parameterizations --\n")
 print (res2.tab)


}
\keyword{models}
\keyword{ts}