`logL.SameMs` <-
function (p, y)
# computes logL over >1 sequence, of model in which all sequences have the 
# same directionality (Mstep), with different step variances
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m, v-1,..v-nseq}
{
  Sm<-0
  nseq<- length(y)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[1],p[i+1]), y[[i]])
   	
  return(Sm) 	
}

