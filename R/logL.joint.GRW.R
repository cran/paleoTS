logL.joint.GRW <-
function (p, x)
# returns logL of GRW model for paleoTS object x
# p is vector of parameters: alpha, ms, vs
{
 # prepare calculations
 anc<- p[1]
 ms<- p[2]
 vs<- p[3]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- vs*outer(x$tt, x$tt, FUN=pmin)
 diag(VV)<- diag(VV) + x$vv/x$nn 
 
 # compute logL based on multivariate normal
 M<- rep(anc, n) + ms*x$tt
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
 
 return(S)		
}
