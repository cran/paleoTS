logL.covTrack <-
function(p, y, z)
# z is covariate; pars = b (slope), evar (variance)
# IMPT: z is of length ns-1; one for each AD transition IN ORDER
{
  b<- p[1]
  evar<- p[2]
  dy <- diff(y$mm)
  dt <- diff(y$tt)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svA <- sv[1:nd]
  svD <- sv[2:(nd + 1)]
  svAD <- svA + svD
  #S <- -0.5 * log(2 * pi * (evar + svAD)) - ((dy - (b*z))^2)/(2 * (evar + svAD))
  S<- dnorm(x=dy, mean=b*z, sd=sqrt(evar+svAD), log=TRUE)
  return(sum(S))	
}
