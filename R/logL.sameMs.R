`logL.sameMs` <-
function (p, y, pool=TRUE)
# computes logL over >1 sequence, of model in which all sequences have the 
# same directionality (Mstep), with different step variances
# y is list of K paleoTS objects, p is array of K+1 parameters {m, v1,..vk}
{
  if (class(y)=="paleoTS")
  	stop("Function logL.sameMs() is only meaningful for multiple sequences.\n")
  	
  Sm<-0
  K<- length(y)
  for (i in 1:K)
   	 Sm<- Sm + logL.RW(p=c(p[1],p[i+1]), y[[i]], pool=pool)
   	
  return(Sm) 	
}

