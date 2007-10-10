`logL.Mult` <-
function (p, y, model=c("GRW", "URW"))
# calculate logL over multiple sequences
#  here, y is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(y)
  for (i in 1:nseq)
    {
     if (model=="URW")
          Smult<- Smult + logL.URW(p,y[[i]])
     else if (model=="GRW")
          Smult<- Smult + logL.GRW (p,y[[i]])
    }
  return (Smult)
}

