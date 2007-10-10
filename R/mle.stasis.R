`mle.Stasis` <-
function (y)
# analytical solution to stasis model
{
  ns<- length(y$mm)
  vp<- pool.var(y)
  th<- mean(y$mm[2:ns])
  om<- var(y$mm[2:ns]) - vp/mean(y$nn)
  
  w<- c(th, om)
  names(w)<- c("theta", "omega")
  return(w)
}

