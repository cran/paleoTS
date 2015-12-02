\name{sim.OU}
\alias{sim.OU}
\alias{ou.M}
\alias{ou.V}

\title{ Simulate evolutionary time-series }
\description{
Generates an evolutionary time-series according to an Orstein-Uhlenbeck (OU) model. }

\usage{
sim.OU(ns = 20, anc = 0, theta = 10, alpha = 0.3, vs = 0.1, vp = 1, 
       nn = rep(20, ns), tt = 0:(ns-1))

ou.M(anc, theta, aa, tt)
ou.V(vs, aa, tt)
}

\arguments{
  \item{ns}{ number of samples in time-series }
  \item{anc}{ ancestral phenotype at the start of the series }
  \item{theta}{ phenotype of the evolutionary optimum }
  \item{alpha}{ strength of the attracting force pulling the population to the optimum }
  \item{vs}{ step variance of the random walk component of change }
  \item{vp}{ within-population trait variance }
  \item{nn}{ vector of the number of individuals in each sample }
  \item{tt}{ vector of sample ages, increases from oldest to youngest }
  \item{aa}{ strength of the attracting force pulling the population to the optimum (same as \code{alpha}) }
}
\details{
See Hansen (1997) for a description of this model in a macroevolutionary context.  This model also arises naturally in 
microevolution as a finite population evolving in the vicinity of an optimum in the adaptive landscape; see Lande (1976) 
and Estes & Arnold (2007). \cr
Functions \code{ou.M} and \code{ou.V} are used internally by \code{sim.OU} to generate the 
means and variances of an OU process.
}

\value{
 A \code{paleoTS} object for \code{sim.OU}.  For \code{ou.M} and \code{ou.V}, a vector of means or variances, respectively, are generated. 
}

\references{ 
Lande, R. 1976. Natural selection and random genetic drift in phenotypic evolution. \emph{Evolution} \bold{ 30}:314-334. \cr
Hansen, T. 1997.  Stabilizing selection and the comparative analysis of adaptation. \emph{Evolution} \bold{ 51}:1341-1351. \cr
Estes, S. & Arnold, S. J. 2007.  Resolving the paradox of stasis: models of stabilizing selection explain evolutionary divergence 
on all timescales. \emph{American Naturalist} \bold{ 169}:227-244. \cr
Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{ 62}:700--710. \cr  
 }
\author{ Gene Hunt }

\seealso{ \code{\link[paleoTS]{sim.GRW}}, \code{\link{opt.joint.OU}} }

\examples{
x1<- sim.OU(ns=100, anc=0, theta=10, alpha=0.2, vs=0.1, vp=0.1, nn=rep(100, times=100), tt=0:99)
plot(x1)
}

\keyword{models}
\keyword{ts}