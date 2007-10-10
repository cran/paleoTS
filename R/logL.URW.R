`logL.URW` <-
function(p,y)
# function to return log-likelihood of step mean and variance (M,V)= p,
# given a paleoTS object
{
  # get parameter values: V is Vstep
  V<- p[1]
    
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  S<- -0.5*log(2*pi*(V*dt+svAD)) - (dy^2)/(2*(V*dt + svAD))
  return(sum(S))
}

