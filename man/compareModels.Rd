\name{compareModels}
\alias{compareModels}

\title{
Compare output from any set of model fits
}

\description{
Function compares a series of objects of \code{paleoTSfit} output from any of the model-fitting functions.
}

\usage{
compareModels(..., silent = FALSE)
}


\arguments{
  \item{\dots}{
a series of one or more \code{paleoTSfit} objects, separated by commas
}
  \item{silent}{
if TRUE, results of model comparison will not be printed to the screen
}
}

\details{
This function offers a convenient way to view the performance of multiple fitted models.  Log-likeihoods, numbers of parameters (K), AICc scores and Akaike weights are returned in a table.
The function checks that all model fits use the same method (parameterization; AD or Joint).  If they do not, an error message results.  Model fits are only comparable if they are based on the same parameterization.
}

\value{
If \code{silent=FALSE}, returns a dataframe with the above information.  If \code{silent=TRUE}, returns a list of the dataframe summarizing model fits along with the maximum-likelihood estimates for the parameters in each model.
}

\author{
Gene Hunt
}



\seealso{
\code{\link{fit3models}}
}

\examples{
	x<- sim.GRW(ns=40)
	m1<- opt.URW(x)
	m2<- opt.GRW(x)
	m3<- fitGpunc(x, ng=2, minb=10)
	compareModels(m1,m2,m3)
}

\keyword{models}
\keyword{ts}
