\name{Stickleback data}
\alias{dorsal.spines}
\alias{pterygiophores}
\alias{pelvic.score}

\docType{data}

\title{ Stickleback data from Bell et al. (2006) }

\description{
Dorsal spine, pterygiophore, and pelvic score data from a fossil stickleback lineage. 
}

\usage{
data(dorsal.spines)
data(pterygiophores)
data(pelvic.score)
	 }

\format{
 \code{paleoTS} objects for time-series of stickleback data.
}

\details{
These data are already in the form of a \code{paleoTS} object, with sample means (\code{mm}), variances (\code{vv})
sample sizes (\code{nn}), and ages in Kyr (\code{tt}).
}

\source{
 Bell, M.A., M.P. Travis and D.M. Blouw 2006. Inferring natural selection in a fossil threespine stickleback. \emph{Paleobiology} \bold{32}:562-577.
}

\examples{
data(dorsal.spines)

# get subset of samples from invading lineage (tt>=4.5 Kyr), only those with nn>=5
ok<- dorsal.spines$tt >= 4.5 & dorsal.spines$nn >=5
ds2<- sub.paleoTS(dorsal.spines, ok=ok)

# convert time scale to generations (500 gen per Kyr)
ds2$tt<- ds2$tt*(1000/2)	
plot(ds2, pool=TRUE)

# compare URW (drift) and OU (adaptive) models
m.urw<- opt.joint.URW(ds2, pool=TRUE)
m.ou<- opt.joint.OU(ds2, pool=TRUE)
aw<- akaike.wts(c(m.urw$AICc, m.ou$AICc))
cat ("Akaike wts: ", round(aw, 3),"\t", "[URW, OU]", "\n")
}

\keyword{datasets}
