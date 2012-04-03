logL.Mult.covTrack <-
function (p, yl, zl)
# y is a _list_ of paleoTS objects, and z is a _list_ of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
    {	Smult<- Smult + logL.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)	
}
