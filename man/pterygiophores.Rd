\name{Stickleback data}
\alias{pterygiophores}
\docType{data}

\title{ Stickleback data from Bell et al. (2006) }

\description{
 Counts of the number of touching pterygiophores from a fossil stickleback lineage. 
}

\usage{data(dorsal.spines)}

\format{
a \code{paleoTS} object for the time-series of pterygiophore data
}

\details{
These data are already in the form of a \code{paleoTS} object, with sample means (\code{mm}), variances (\code{vv})
sample sizes (\code{nn}), and ages in Kyr (\code{tt}).
}

\source{
 Bell, M.A., M.P. Travis and D.M. Blouw 2006. Inferring natural selection in a fossil threespine stickleback. \emph{Paleobiology} \bold{32}:562-577.
}

\examples{
data(pterygiophores)

# get subset of samples from invading lineage (tt>=4.5 Kyr), only those with nn>=5
ok<- pterygiophores$tt >= 4.5 & pterygiophores$nn >=5
pt2<- sub.paleoTS(pterygiophores, ok=ok)

# convert time scale to generations (500 gen per Kyr)
pt2$tt<- pt2$tt*(1000/2)	
plot(pt2, pool=TRUE)

# compare URW (drift) and OU (adaptive) models
m.urw<- opt.joint.URW(pt2, pool=TRUE)
m.ou<- opt.joint.OU(pt2, pool=TRUE)
aw<- akaike.wts(c(m.urw$AICc, m.ou$AICc))
cat ("Akaike wts: ", round(aw, 3),"\t", "[URW, OU]", "\n")

}

\keyword{datasets}
