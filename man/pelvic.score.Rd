\name{Stickleback data}
\alias{pelvic.score}

\docType{data}

\title{ Stickleback data from Bell et al. (2006) }

\description{
 Data of pelvic score data from a fossil stickleback lineage. 
}

\usage{data(dorsal.spines)}

\format{
a \code{paleoTS} object for the time-series of pelvic scores
}

\details{
These data are already in the form of a \code{paleoTS} object, with sample means (\code{mm}), variances (\code{vv})
sample sizes (\code{nn}), and ages in Kyr (\code{tt}).
}

\source{
 Bell, M.A., M.P. Travis and D.M. Blouw 2006. Inferring natural selection in a fossil threespine stickleback. \emph{Paleobiology} \bold{32}:562-577.
}

\examples{
data(pelvic.score)

# get subset of samples from invading lineage (tt>=4.5 Kyr), only those with nn>=5
ok<- pelvic.score$tt >= 4.5 & pelvic.score$nn >=5
ps2<- sub.paleoTS(pelvic.score, ok=ok)

# convert time scale to generations (500 gen per Kyr)
ps2$tt<- ps2$tt*(1000/2)	
plot(ps2, pool=TRUE)
}

\keyword{datasets}
