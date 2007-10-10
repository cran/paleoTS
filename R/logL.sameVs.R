`logL.SameVs` <-
function (p, y)
# computes logL over >1 sequence, of model in which all sequences have the 
# same step variance, with different steo means
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m-1,..m-nseq, vs}
{
  Sm<-0
  nseq<- length(y)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[i],p[nseq+1]), y[[i]])
   	
  return(Sm) 	
}

