\name{paleoTS-package}
\alias{paleoTS-package}
\alias{paleoTS}
\docType{package}
\title{
Analysis of evoltuionary time-series 
}
\description{
This package facilitates analysis of paleontological sequences of trait values from an evolving lineage.  Functions are provided to fit, using maximum likelihood, evolutionary models including unbiased random walks, directional evolution, stasis, Ornstein-Uhlenbeck, punctuated change, and evolutionary models in which traits track some measured covariate.
}
\details{
\tabular{ll}{
Package: \tab paleoTS\cr
Type: \tab Package\cr
Version: \tab 0.4-3\cr
Date: \tab 2011-12-22\cr
License: \tab GPL 2\cr
}

}
\author{ Gene Hunt}

\references{
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{ 32}:578--601.
Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{ 62}:700-710. \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology}\bold{ 34}:360--377.  \cr
Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based assessment of modes, rates, punctuations and process..  \emph{In} R.K. Bambach and P.H. Kelley, eds.  From Evolution to Geobiology: Research Questions Driving Paleontology at the Start of a New Century:578--601. \cr
Hunt, G., S. Wicaksono, J. E. Brown, and G. K. Macleod. 2010. Climate-driven body size trends in the ostracod fauna of the deep Indian Ocean. \emph{Palaeontology} 53(6):1255-1268. \cr

}

\keyword{ package }
\examples{
 x <- sim.GRW(ns=30, ms=0, vs=0.3)	# simulate unbiased random walk
 fit3models(x)		# compare fits of directional, random walk, and stasis models
}
