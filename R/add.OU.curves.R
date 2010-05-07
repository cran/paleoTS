add.OU.curves <-
function(w, x, what=c("lines", "polygon"), tt.offset=0, ...)
# add expected OU trajectory to a plot
# w is output from opt, x is paleoTS
{
 tt2<- x$tt-min(x$tt)
 ee<- ou.M(w$par[1], w$par[3], w$par[4], tt2)
 vv<- ou.V(w$par[2], w$par[4], tt2)
 se<- sqrt(vv)
 
 what<- match.arg(what)
 
 if (what=="lines")
 {
  lines(x$tt+tt.offset, ee, ...)
  lines(x$tt+tt.offset, ee+2*se, ...)
  lines(x$tt+tt.offset, ee-2*se, ...)
 }
 if(what=="polygon"){	
  polygon(c(x$tt+tt.offset, rev(x$tt+tt.offset)), c(ee+2*se,rev(ee-2*se)), ...)
  lines (x$tt+tt.offset, ee, col="black", lwd=1)
  }
}

