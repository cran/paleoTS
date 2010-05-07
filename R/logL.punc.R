logL.punc <-
function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S+ logL.Stasis(p=c(th[i], om[i]), xl[[i]])
  
  return (S)	
}

