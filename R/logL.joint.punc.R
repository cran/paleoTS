logL.joint.punc <-
function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  M<- th[gg]  # vector of MVN means
  VV<- diag(om[gg] + y$vv/y$nn)  # vcv matrix
  S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  
  return (S)	
}
