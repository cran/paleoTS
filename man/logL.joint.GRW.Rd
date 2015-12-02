\name{logL.joint.GRW}
\alias{logL.joint.GRW}
\alias{logL.joint.URW}
\alias{logL.joint.Stasis}
\alias{logL.joint.StrictStasis}
\alias{logL.joint.OU}

\title{ Log-likelihoods for evolutionary models (joint parameterization) }

\description{
Returns log-likelihood for general random walk (\code{logL.joint.GRW}), unbiased random walk (\code{logL.joint.URW}), stasis (\code{logL.joint.Stasis}) and OU (\code{logL.joint.OU}) models for the \code{Joint} parameterization.
}

\usage{
logL.joint.GRW(p, y)
logL.joint.URW(p, y)
logL.joint.Stasis(p, y)
logL.joint.StrictStasis(p, y)
logL.joint.OU(p, y)
}

\arguments{
  \item{p}{ a vector of parameters }
  \item{y}{ a \code{paleoTS} object }
}

\details{
 For the general random walk, \code{p = c(anc, mstep, vstep)}; for an unbiased random walk, \code{p = c(anc, vstep)}; for the stasis model, \code{p = c(theta, omega)}, and for the OU model \code{p = c(anc, vstep, theta, alpha)}.  In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values. }

\value{
  The log-likelihood of the parameter estimates (\code{p}), given the data (\code{x}).
  }
  
\references{ 
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601. \cr
	Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{ 62}:700-710. \cr   
	Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based assessment of modes, rates, punctuations and process..  \emph{In} R.K. Bambach and P.H. Kelley, eds.  From Evolution to Geobiology: Research Questions Driving Paleontology at the Start of a New Century:578--601. \cr
	}
\author{ Gene Hunt }

\note{
	The user will not generally need to call these functions directly; they are used by the optimization functions.
}

\section{Warning}{
Because these functions parameterize the models differently, their log-likelihoods are not comparable to those that do not use the joint parameterization.
}

\seealso{ \code{\link{opt.joint.GRW}}, \code{\link{logL.GRW}} }

\examples{
x<- sim.GRW(ns=20, ms=0, vs=1)
L1<- logL.joint.GRW(p=c(0,0,1), x)	# actual parameters
L2<- logL.joint.GRW(p=c(0,10,1), x)	# should be a bad guess
cat(L1, L2, "\n")
}

\keyword{models}
\keyword{ts}
