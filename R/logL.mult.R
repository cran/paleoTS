`logL.mult` <-
function (p, y, model=c("RW", "RWu"), pool=TRUE)
# calculate logL over multiple sequences
#  here, y is a list of paleoTS sequences
#  variances are pooled (if desired) separately for each sequence
{
  if (class(y)=="paleoTS")
     stop("Use logL.RW() or logL.RWu for single paleoTS sequences\n")

  model<- match.arg(model)
  Smult<-0
  nseq<- length(y)
  for (i in 1:nseq)
    {
     if (model=="RWu")
          Smult<- Smult + logL.RWu(p,y[[i]], pool=pool)
     else if (model=="RW")
          Smult<- Smult + logL.RW (p,y[[i]], pool=pool)
    }
  return (Smult)
}

