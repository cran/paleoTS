logL.joint.punc.omega <-
function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  M<- th[gg]  # vector of MVN means
  omv<- rep(om, max(gg))
  VV<- diag(omv[gg] + y$vv/y$nn)  # vcv matrix
  S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  
  return (S)	  	
}
