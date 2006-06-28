"plot.paleoTS" <-
function (x, nse=1, yl=NULL, pch=19, true.means=FALSE, start.time=NULL, add=FALSE, ...)
# plots paleoTS object, with nse*se error bars
{
  se<- sqrt(x$vv/x$nn)
  lci<- x$mm-(nse*se)
  uci<- x$mm+(nse*se)
  if (is.null(yl))
     yl<-range(c(lci, uci))

  if (!is.null(start.time))
     x$tt<- x$tt-start.time
     
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
    #mtext(x$label, side=3, col="grey", font=3, cex=0.6)
   }
  segments (x$tt, lci, x$tt, uci, lty=1, ...)
  
  if (true.means)
  	lines (x$tt, x$MM, col="blue", lty=2)
}

