logL.joint.URW <-
function(p,x)
# returns logL of URW model for paleoTS object x
# p is vector of parameters: alpha, vs
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- vs*outer(x$tt, x$tt, FUN=pmin)
 diag(VV)<- diag(VV) + x$vv/x$nn 	# add sampling variance
 detV<- det(VV)
 invV<- solve(VV)
 
 # compute logL based on multivariate normal
 M<- rep(anc, n)
 S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)

 return(S)
}

