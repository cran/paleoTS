`read.paleoTS` <-
function (file=NULL, hh=FALSE, oldest="first", ...)
# read in paleoTS data from a file
# samples should be listed from oldest (first) to youngest (last)
# if time is in ages BP (oldest with largest age), convert to forward moving time
#    by tt.new<- max(tt)-tt.old
{
  if (is.null(file))
     {
      ff<- file.choose()
      x<-read.table(ff, header=hh, ...)
      lab1<- ff
     }

  else
     {
       x<-read.table(file=file, header=hh, ...)
       lab1<-paste(getwd(), file)
     }
      
  # change to proper paleoTS object
  xr<- as.paleoTS(mm=x[,2], vv=x[,3], nn=x[,1], tt=x[,4], label=lab1)


  if (oldest=="last")
    {
      # reverse order of samples	
      oo<- length(xr$mm):1
      xr$mm <- xr$mm[oo]
      xr$vv <- xr$vv[oo]
      xr$nn <- xr$nn[oo]
      xr$tt <- xr$tt[oo]	
    }
  
  if (xr$tt[1]>xr$tt[2])	# ages are in Ka/Ma etc, so numbers decrease
  	{   xr$start.age<- xr$tt[1]
  		xr$tt<- max(xr$tt)-xr$tt }			# change ages so oldest is 0
  		
  return (xr)
}

