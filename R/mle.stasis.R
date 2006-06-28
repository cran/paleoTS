`mle.stasis` <-
function (x)
# analytical solution to stasis model
{
  ns<- length(x$mm)
  vp<- pool.var(x)
  th<- mean(x$mm[2:ns])
  om<- var(x$mm[2:ns]) - vp/mean(x$nn)
  
  w<- c(th,om)
  names(w)<- c("theta.hat","omega.hat")	
  return(w)
}

