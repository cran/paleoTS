`split4punc` <-
function (y, gg, overlap=TRUE)
# divides a paleoTS object (y) into a several paleoTS objects, according to vector 'gg'
# gg is a vectors of 1,2,3.. indicating groupings
# overlap=TRUE means that the adjacent samples are included 
{
  yl<- list()
  ng<- max(gg)
  for (i in 1:ng)
   {
   	 ok<- gg==i
   	 if (i<ng & overlap==TRUE)
   	 	ok[min(which(gg==i+1))]<- TRUE		# tack on next sample too
   	
   	 yl[[i]]<- as.paleoTS(y$mm[ok],y$vv[ok],y$nn[ok],y$tt[ok],y$MM[ok])
   }
  return (yl)	
}

