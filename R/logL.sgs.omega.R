logL.sgs.omega <-
function(p, y, gg, model="GRW")
## log-likelihood of sgs model, omega shared over stasis segments
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1]
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5]	}
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4]	}
  
  l1<- logL.Stasis(p=c(th[1], om), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om), yl[[3]])
  
  logl<- l1+l2+l3
  return(logl)  		
}
