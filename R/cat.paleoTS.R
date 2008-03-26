`cat.paleoTS` <-
function (y)
# concatenates multiple paleoTS objects, with y a list of paleoTS objects
{
 x<- y[[1]]
 for (i in 2:length(y))
   {
   	x$mm<- append(x$mm, y[[i]]$mm)
   	x$vv<- append(x$vv, y[[i]]$vv)
   	x$tt<- append(x$tt, y[[i]]$tt)
   	x$MM<- append(x$MM, y[[i]]$MM)
   	x$nn<- append(x$nn, y[[i]]$nn)
   }
  
  return (x)   	
}

