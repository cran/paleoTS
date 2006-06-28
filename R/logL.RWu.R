"logL.RWu" <-
function(p,y, pool=TRUE)
# function to return log-likelihood of step mean and variance (V)=p
# constrained to be unbiased RW (M=0)
# given a paleoTS object
{
   # calculate mean and age differences
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  vp<- pool.var(y)

  V<- p
  S<-0
  for (i in 1:nd)
    {
      if (pool)
      	Sadd<- -0.5*log(2*pi*(V*dt[i]+(vp/y$nn[i])+(vp/y$nn[i+1]))) - (dy[i]^2)/(2*(V*dt[i]+(vp/y$nn[i])+(vp/y$nn[i+1])))
      else
      	Sadd<- -0.5*log(2*pi*(V*dt[i]+(y$vv[i]/y$nn[i])+(y$vv[i+1]/y$nn[i+1]))) - (dy[i]^2)/(2*(V*dt[i]+(y$vv[i]/y$nn[i])+(y$vv[i+1]/y$nn[i+1])))
      S<- S + Sadd
    }

  return (S)
}

