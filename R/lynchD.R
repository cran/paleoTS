`lynchD` <-
function (y, gen.per.t=1e6, pool=TRUE, ...)
# compute Lynch's rate metric, Delta
# gen.per.t is the number of generations per unit tt
# ... are further arguments to opt.URW()
{
  w<- opt.URW(y, pool=TRUE, ...)
  vp<- pool.var(y)
  D<- (unname(w$par)/vp)/gen.per.t
  
  drift.min<- 5e-5
  drift.max<- 5e-3
  
  if (D < drift.min)	res<- "Slower than drift range"
  else if (D > drift.max)	res<- "Faster than drift range"
  else					res<- "Within range of drift"
  
  w<- list(D=D, pooled.var=vp, gen.per.t=gen.per.t, vstep=w$par, drift.range=c(drift.min, drift.max), result=res) 
  return(w)	
}

