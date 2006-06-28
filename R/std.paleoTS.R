`std.paleoTS` <-
function (y, zero="start")
# returns paleoTS, converted in phenotypic SD units
# mm -> (mm - mean(mm) )/sqrt(vp); vv -> vv/vp
# optionally set starting mean value to zero
{
 vp<- pool.var(y)
 sp<- sqrt(vp)
 
 ys <- y 
 ys$mm<- (y$mm - mean(y$mm)) /sp
 ys$vv<- y$vv/vp
 
 if (zero=="start")
 	ys$mm <- ys$mm - ys$mm[1]
 
 return(ys)	
}

