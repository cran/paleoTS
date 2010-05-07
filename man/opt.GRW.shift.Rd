\name{opt.GRW.shift}
\alias{opt.GRW.shift}

\title{ Functions for random walks with shifting parameters }

\description{
Functions to simulate and to infer a model with random walk dynamics, with parameter values that 
shift at one or more points in the sequence.
}

\usage{
opt.GRW.shift(y, ng = 2, minb = 5, model = 1, pool = TRUE, silent = FALSE)
}

\arguments{
  \item{y}{ a \code{paleoTS} object }
  \item{ng}{ the number of different segments in the sequence }
  \item{minb}{ the minimum number of samples to consider as a segment }
  \item{model}{ options for variants of random walk to fit (see \code{Details}). }
  \item{pool}{ if TRUE, pool phenotypic variances across samples }
  \item{silent}{ if TRUE, do not print information on fitting to screen }
}
\details{
This model divides an evolutionary sequence into two or more non-overlapping parts 
called segments, and then fits some version of the random walk model to each segment separately. \cr
The \code{model} argument has four options: \cr
  \code{model=1}, general random walk model, step variance shared across segments \cr
  \code{model=2}, general random walk model, step mean shared across segments \cr
  \code{model=3}, unbiased random walk, separate step variance for each segment
  \code{model=4}, general random walk, speparate step mean and variance for each segment
}
\value{
  A list including:
  \item{value }{ the log-likelihood of the optimal solution}
  \item{par }{ parameter estimates}
  \item{K}{ number of parameters in the model}
  \item{n}{ the number of observations, equal to the number of evoltuionary transistions}
  \item{shift.start}{ the index of the initial samples of each segment}
  \item{AIC}{ Akaike information criterion}
  \item{AICc}{ modified Akaike information criterion}
  \item{BIC}{ Bayes (or Schwarz) information criterion}
  \item{shift.start}{ index of each sample that starts a new segment}
  \item{all.logl}{ log-likelihoods for all tested partitions of the series into segments}
  \item{GG}{ matrix of indices of initial samples of each tested segment configuration; each column of \code{GG} corresponds to the elements of \code{all.logl}}
}
\references{ 
Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond.  \emph{Paleobiology} \bold{32}:578--601.  \cr
Hunt, G. 2008. Gradual or pulsed evolution: when should punctuational explanations be preferred?  \emph{Paleobiology} \bold{34}:360--377.  \cr
}
\author{ Gene Hunt}

\seealso{ \code{\link{sim.GRW.shift}}, \code{\link{opt.GRW}}, \code{\link{opt.RW.Mult}} }
\examples{
 x<- sim.GRW.shift(ns=c(20,20), ms=c(0,1), vs=c(0.2, 0.2))
 plot(x)
 w.shift<- opt.GRW.shift(x, ng=2, model=1)
 print (w.shift$par)
 print (w.shift$shift.start)  # estimated first sample in second segment
}

\keyword{models}
\keyword{ts}