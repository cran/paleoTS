`plot.paleoTS` <-
function (x, nse=1, yl=NULL, pch=NULL, pool=TRUE, true.means=FALSE, add=FALSE, ...)
# plots paleoTS object, with nse*se error bars
{
  if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
  se<- sqrt(x$vv/x$nn)
  lci<- x$mm-(nse*se)
  uci<- x$mm+(nse*se)
  if (is.null(yl))
     yl<-range(c(lci, uci))

  if (!is.null(x$start.age))
     x$tt<- x$tt-x$start.age
     
  if (is.null(pch))
  	pch=19
  
  if (add)
   {
     lines (x$tt, x$mm, ...)
     points (x$tt, x$mm, ...)
   }
  else
   {
    plot (x$tt, x$mm, type="o", pch=pch, xlab="Time", ylab="Trait Mean", ylim=yl, ...)
   }
  segments (x$tt, lci, x$tt, uci, lty=1, ...)
  
  if (true.means)
  	lines (x$tt, x$MM, col="blue", lty=2)
}

