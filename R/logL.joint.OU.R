logL.joint.OU <-
function(p,x)
# returns logL of OU model for paleoTS object x
# 
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 theta<- p[3]
 aa<- p[4]
 n<- length(x$mm)
 
 # compute covariance matrix
 ff<- function (a,b) abs(a-b)
 VV<- outer(x$tt, x$tt, FUN=ff)
 VV<- exp(-aa*VV)
 VVd<- ou.V(vs,aa,x$tt)
 VV2<- outer(VVd,VVd,pmin)
 VV<- VV*VV2
 diag(VV)<- VVd + x$vv/x$nn 	# add sampling variance
 detV<- det(VV)
 invV<- solve(VV)
 
 # compute logL based on multivariate normal
 M<- ou.M(anc, theta, aa, x$tt)
 S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 return(S)		
}

