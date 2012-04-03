logL.punc.omega <-
function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S + logL.Stasis(p=c(th[i], om), xl[[i]])
  
  return (S)	  	
}
