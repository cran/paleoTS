`test.var.het` <-
function (y, method="Bartlett")
# test for variance heterogeneity among samples in a paleoTS object
{

 vp<- pool.var(y)
 NN<- sum(y$nn)
 k<- length(y$vv)
 top<- (NN-k)*log(vp)-(sum((y$nn-1)*log(y$vv)))
 bot<- 1+ (1/(3*(k-1)))*( (sum(1/(y$nn-1))) - (1/(NN-k)))
 TT<- top/bot
 p.val<- pchisq(TT, df=k-1, lower.tail=FALSE)

 w<-list(stat=TT, p.value=p.val, df=k-1)
 return (w)
}

