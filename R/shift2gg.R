shift2gg <-
function (ss, ns)
# ss is vector of shift points, ns is # samples
{
  z<- c(0,ss,ns+1)
  cc<- cut(1:ns, breaks=z, right=FALSE)
  gg<- as.numeric(cc)
  return(gg)	
}
