`ln.paleoTS` <-
function (y)
# returns paleoTS, with data approx ln-transformed
# mean(ln[y])= ln(mean[y]); var(ln[y])=(sd[y]/mean[y])^2
{
 logx<- y
 logx$mm<- log(y$mm)
 logx$vv<- (sqrt(y$vv)/y$mm)^2
 
 return (logx)
}

