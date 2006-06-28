`logL.stasis` <-
function(p, y, pool=TRUE)
## logL of stasis, same white noise as stasis1, but based on differences in mean
## to be comparable to the other models. 
## Under this model, dx is not random draw, but depends on the value at the start
## dx ~ N(M-xi, Ve + Vs)	where:
##		xi is the phenotypic value at the start of the step
##		M is the central peak (free parameter)
##		Ve is the 'evolutionary' variance (free parameter)
##		Vs is the sampling variance of dx		
## time does not enter logL; it is non-Markovian
{
 M<-p[1]
 V<-p[2]
 dy<- diff(y$mm)
 nd<- length(dy)

 # compute sampling variances for each sample
 if (pool)
 	{
 	 vp<- pool.var(y)
 	 sv<- rep(vp/y$nn, nd+1)			
 	}
 else	sv<- y$vv/y$nn
 
 S<-0
 for (i in 1:nd)
 	{
	  #svi<- sv[i] + sv[i+1]	# don't need sv[i], bc conditional on observed xi
	  svi<- sv[i+1]
	  Sadd<- -0.5*log(2*pi*(V+svi)) - (dy[i]-(M-y$mm[i]))^2/(2*(V+svi)) 
	  S<- S + Sadd			
 	}	
  return (S)	

}

