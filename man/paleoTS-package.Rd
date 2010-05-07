\name{paleoTS-package}
\alias{paleoTS-package}
\alias{paleoTS}
\docType{package}
\title{
Analysis of evoltuionary time-series 
}
\description{
This package includes a variety of functions useful for analyzing evoltuionary time-series, 
particularly for paleontological data.  The functions are useful for analyzing sequences in terms 
of statistical models, with three models implemented: directional evolution (modeled as a general random walk),
unbiased random walks, and stasis.  These analyses are outlined in Hunt (2006).  Two different parameterizations 
can be used to fit these models, the default (based on ancestor-descendant differences), and a joint method, which 
uses all sample means jointly to fit the model.  The joint parameterization can also be used to fit an Ornstein-Uhlenbeck model, 
which is intermediate between a random walk and stasis (see Hunt et al., 2008). \cr
Starting with v0.3-1, functions are included that allow fitting of models of punctuated evolution, and other models in which 
evolutionary dynamics change throughout an observed evolutionary sequence (see Hunt 2008).
}
\details{
\tabular{ll}{
Package: \tab paleoTS\cr
Type: \tab Package\cr
Version: \tab 0.3-1\cr
Date: \tab 2008-01-16\cr
License: \tab GPL 2\cr
}

}
\author{ Gene Hunt}

\references{
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601.
Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{62}:700-710. \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology}\bold{34}:360--377.  \cr
Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based assessment of modes, rates, punctuations and process..  \emph{In} R.K. Bambach and P.H. Kelley, eds.  From Evolution to Geobiology: Research Questions Driving Paleontology at the Start of a New Century:578--601. \cr
}

\keyword{ package }
\examples{
 x <- sim.GRW(ns=30, ms=0, vs=0.3)	# simulate unbiased random walk
 fit3models(x)		# compare fits of directional, random walk, and stasis models
}
