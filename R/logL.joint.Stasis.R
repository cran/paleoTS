`logL.joint.Stasis` <-
function (p,x)
# returns logL of Stasis model for paleoTS object x
# p is vector of parameters: theta, omega
{
 # prepare calculations
 theta<- p[1]
 omega<- p[2]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- diag(omega + x$vv/x$nn)  # omega + sampling variance
 detV<- det(VV)
 invV<- solve(VV)
 
 # compute logL based on multivariate normal
 M<- rep(theta, n)
 S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 return(S)	 	
}

