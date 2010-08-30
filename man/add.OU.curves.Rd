\name{add.OU.curves}
\alias{add.OU.curves}

\title{ Adds curves to an existing plot indicating best-fit OU model }
\description{
This function takes the results of fitting the Orstein-Uhlenbeck (OU) model and plots the expected
 value (mean) and 95 percent probability interval of the best-fitting OU parameters.
}
  
\usage{
add.OU.curves(w, x, what = c("lines", "polygon"), tt.offset = 0, ...)
}

\arguments{
  \item{w}{ a paleoTSfit object output from the function \code{opt.joint.OU} }
  \item{x}{ a \code{paleoTS} pbject }
  \item{what}{ whether to draw \code{lines} or \code{polygon} for the model fit }
  \item{tt.offset}{ difference between starting age in \code{w} and \code{x} }
  \item{\dots}{ additional arguments to graphing functions }
}
\details{
Note that \code{tt.offset} is only used in the relatively rare instances in which the model fit 
is not applied to the entire sequence.
}
\value{
None, the function is evoked to add graphics to an existing plot.
}

\author{ Gene Hunt }

\seealso{ \code{\link[paleoTS]{plot.paleoTS}}, \code{\link{opt.joint.OU}} }

\examples{
data(dorsal.spines)
# get subset of samples from invading lineage (tt>=4.5 Kyr), only those with nn>=5
ok<- dorsal.spines$tt >= 4.5 & dorsal.spines$nn >=5
ds2<- sub.paleoTS(dorsal.spines, ok=ok)

# convert time scale to generations (500 gen per Kyr)
ds2$tt<- ds2$tt*(1000/2)	
plot(ds2, pool=TRUE)

# fit OU model
m.ou<- opt.joint.OU(ds2, pool=TRUE)
add.OU.curves(m.ou, ds2, what="lines", col="darkgrey")

}

\keyword{models}
\keyword{ts}
