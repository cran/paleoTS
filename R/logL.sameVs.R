`logL.sameVs` <-
function (p, y, pool=TRUE)
# computes logL over >1 sequence, of model in which all sequences have the 
# same step variance, with different steo means
# y is list of K paleoTS objects, p is array of K+1 parameters {m1,..mk, vs}
{
  if (class(y)=="paleoTS")
  	stop("Function logL.sameVs() is only meaningful for multiple sequences.\n")
  	
  Sm<-0
  K<- length(y)
  for (i in 1:K)
   	 Sm<- Sm + logL.RW(p=c(p[i],p[K+1]), y[[i]], pool=pool)
   	
  return(Sm) 	
}

