`as.paleoTS` <-
function (mm, vv, nn, tt, MM=NULL, genpars=NULL, label="")
# converts vectors of observed mean, variances, N's, ages (time steps),
# true means (for simulated series), and label to an object of class 'paleoTS'
{
  res<- list(mm=mm,vv=vv,nn=nn,tt=tt,MM=MM,genpars=genpars,label=label)
  class(res)<- "paleoTS"
  return (res)
}

