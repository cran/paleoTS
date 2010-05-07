\name{opt.joint.GRW}
\alias{opt.joint.GRW}
\alias{opt.joint.URW}
\alias{opt.joint.Stasis}
\alias{opt.joint.OU}

\title{ Optimize evolutionary models (joint parameterization) }

\description{
Functions to find maximum likelihood solutions to general random walk (\code{opt.joint.GRW}), unbiased random walk (\code{opt.joint.URW}), stasis (\code{opt.joint.Stasis}) and OU models (\code{opt.joint.OU}).
}

\usage{
opt.joint.GRW(x, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.joint.URW(x, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.joint.Stasis(x, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
opt.joint.OU(x, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{paleoTS} object }
  \item{pool}{ logical indicating whether to pool variances across samples }
  \item{cl}{ control list, passed to function \code{optim} }
  \item{meth}{ optimization method, passed to function \code{optim} }
  \item{hess}{ logical, indicating whether to calculate standard errors from the Hessian matrix }
}
\details{
  These functions numerically search a log-likelihood surface for its optimum--they are a convenient wrapper to \code{optim}.  
  Arguments \code{meth}, \code{cl}, and \code{hess} are passed to \code{optim}; see the help for that function for details.  
  These are included to allow sophisticated users greater control over the optimization; the defaults seem to work well for most, 
  but not all sequences. For \code{meth="L-BFGS-B"}, some parameters are constrained to be non-negative, which is useful parameters 
  which cannot truly be negative, such as \code{vstep} (random walk) and \code{omega} (stasis model). \cr
  
  Initial estimates to start the optimization come in part from analytical solutions based on assuming equal sampling error across 
  samples and evenly spaced samples in time (functions \code{mle.GRW}, \code{mle.URW} and \code{mle.Stasis}). 

}
\value{
  \item{par }{ parameter estimates}
  \item{value }{ the log-likelihood of the optimal solution}
  \item{counts}{ returned by \code{optim} }
  \item{convergence}{ returned by \code{optim} }
  \item{message}{ returned by \code{optim} }
  \item{p0}{ initial guess for parameter values at start of optimization}
  \item{K}{ number of parameters in the model}
  \item{n}{ the number of observations, equal to the number of samples}
  \item{AIC}{ Akaike information criterion}
  \item{AICc}{ bias-corrected Akaike information criterion}
  \item{BIC}{ Bayes (or Schwarz) information criterion}
  \item{se}{ standard errors for parameter estimates, computed from the curvature of the log-likelihood surface (only if \code{hess = TRUE})}
  \item{\dots}{ other output from call to \code{optim}}
}

\references{ 
	Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601. \cr
	Hunt, G., M. Bell & M. Travis. 2008.  Evolution towards a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} \bold{62}:700-710. \cr   
	Hunt, G. 2008. Evolutionary patterns within fossil lineages: model-based assessment of modes, rates, punctuations and process..  \emph{In} R.K. Bambach and P.H. Kelley, eds.  From Evolution to Geobiology: Research Questions Driving Paleontology at the Start of a New Century:578--601. \cr
	}
\author{ Gene Hunt }

\note{
These optimizations are performed using a parameterization of the GRW, URW and Stasis models that differs from the functions \code{opt.GRW},
\code{opt.URW}, and \code{opt.Stasis}.  Those functions models are fit from the \emph{differences} between adjacent samples, 
removing the autocorrelation in the time-series.  In the joint parameterization implemented in functions \code{opt.alt.GRW}, \code{opt.alt.URW}, 
\code{opt.alt.Stasis} and \code{opt.alt.OU}, models are fit using the actual sample values, with the autocorrelation among samples accounted for in 
the log-likelihood function.  For each model, the joint distribution of sample means is multivariate normal, with means 
and variance-covariances determined by evolutionary parameters and sampling errors.  Note that the Orstein-Uhlenbeck model is at present only implemented 
using the joint parameterization. \cr

For details on this parameterization of the models, see the paper Hunt, et al (2008).  In general, 
the two different parameterizations yield similar results as to the relative support for different models.  In my experience, the two approaches tend to differ appreciably 
only with relatively long sequences that have small differences between consecutive samples.  In these circumstances, the alternative parameterization is better able to 
distinguish true evolutionary patterns from sampling error, and is less prone to falsely favoring the Stasis model. \cr

The general random walk (\code{GRW}) and unbiased random walk (\code{URW}) models require an additional parameter that specifies the phenotype at the start of 
the sequence.  This parameter is called \code{anc} in the vector of returned parameter estimates, and in general is rather similar
to the trait value of the first sample in the sequence.  This extra parameter is not necessary for the stasis model because the initial value does not 
figure into the likelihood calculations (which are comletely determined by the position of optimum and the variance around this optimum).
}

\section{Warning}{
Measures of model fit (log-likelihoods, AIC scores, etc) are not comparable between the two parameterizations.
}

\seealso{ \code{\link{logL.joint.GRW}}, \code{\link[paleoTS]{opt.GRW}} }

\examples{
 x<- sim.GRW(ns=30, ms=1, vs=1)
 plot(x)
 
 # easier to use  fit3models.joint() 
 m.urw<- opt.joint.URW(x)
 m.grw<- opt.joint.GRW(x)
 m.sta<- opt.joint.Stasis(x)
 
 cat(m.urw$AICc, m.grw$AICc, m.sta$AICc, "\n")	# print AICc scores

}
\keyword{models}
\keyword{ts}