`akaike.wts` <-
function(aa)
## aa is a vector of AIC or AICc values
{
 ma<- min(aa)	
 delt<- aa - ma		# delta values
 denom<- sum(exp(-delt/2))
 
 ww<- exp(-delt/2) / denom
 names(ww)<- names(aa)
 return(ww)	
}

