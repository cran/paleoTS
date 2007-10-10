`logL.Stasis` <-
function(p, y)
## logL of stasis model
{
 # get parameter estimates
 M<-p[1]	# M is theta
 V<-p[2]	# V is omega
 dy<- diff(y$mm)
 nd<- length(dy)
 
 sv<- y$vv/y$nn
 svD<- sv[2:(nd+1)]  # only need sampling varian descendant
 anc<- y$mm[1:nd]

 S<- -0.5*log(2*pi*(V+svD)) - ((dy-(M-anc))^2)/(2*(V+svD))
 return(sum(S))
}

