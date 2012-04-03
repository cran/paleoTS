LRI <-
function(x, gen.per.t=1e6, draw=TRUE)
{
  n<- length(x$mm)
  x$tt<- x$tt*gen.per.t  # convert to generational timescale (e.g., 1e6 when time is measured in Myr and generation time = 1 yr)
  sp<- sqrt(pool.var(x))
  dy<- outer(x$mm, x$mm, FUN='-')  # outer() much faster than looping
  dt<- outer(x$tt, x$tt, FUN='-')
  dy<- abs(dy)
  dt<- abs(dt)
   
  dy<- dy/sp
  ut<- upper.tri(dy)
  logI<- log10(dt[ut])
  logR<- log10(dy[ut]) - logI
   
  # subset only non-zero rates
  ok<- is.finite(logR)
  if (sum(ok) != length(logI))	warning("Some zero rates were ignored.")
  
  # function to fit line in LRI plot robustly, using min abs deviations as per Gingerich 1993
  opt.lad<- function (x, y)
	{
	  ok<- is.finite(x) & is.finite(y)
	  xok<- x[ok]
	  yok<- y[ok]
	  w.ls<- lm(yok ~ xok)  # use LS coef as starting point in optimization
	
	  # function minimized for least abs deviation
	  fad<- function(p, x, y)  sum( abs(y - p[2]*x - p[1]) )
	  w.lad<- optim(w.ls$coef, fn=fad, x=xok, y=yok)
	  return(w.lad$par)	
	}
    
  # call robust regression function
  lf<- opt.lad(logI[ok], logR[ok])
  lf[3]<- 10^lf[1]   # generational rate
  names(lf)<- c("Intercept", "slope", "GenerationalRate")
  
  # do LRI plot, if desired
  if (draw)
   {
   	xl<- c(0, max(logI[ok]))
   	yl<- c(min(logR[ok]), max(lf[1], max(logR[ok])))
   	plot(logI[ok], logR[ok], xlim=xl, ylim=yl, xlab='log10 Interval [generations]', ylab='log10 Rate', cex=0.6)
   	abline(lf[1], lf[2], col='black', lwd=2)
   	title("LRI plot")
   	mtext(paste('data label: ', x$lab), cex=0.6, font=3, col='darkgrey')
   	restext<- paste('Slope = ', round(lf[2],3), '\n', 'Intercept = ', round(lf[1],3), '\n', 'Generational Rate = ', round(lf[3],5), '\n', sep='')
   	text(0, min(logR[ok]), restext, adj=c(0,0), cex=0.7, font=2)
   }
  
  #w<- list(b0=b0, b1=b1, h0=h0, logR=logR, logI=logI, dy=dy)  
  return(lf)
}
