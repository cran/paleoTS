`sub.paleoTS` <-
function (y, ok=NULL, k=0.1)
# subsample a paleoTS, either from steps given by T/F vector 'ok'
#  proportion 'k' of samples, chosen randomly
{
 ys<- y
 ns<- length(y$mm)
 take<- array(FALSE, dim=ns)
 if (!is.null(ok) )
   take<- ok
 else
   take[sample(1:ns, size=round(k*ns))]<- TRUE

 ys$mm<- y$mm[take]
 ys$vv<- y$vv[take]
 ys$nn<- y$nn[take]
 ys$tt<- y$tt[take]
 ys$MM<- y$MM[take]
 ys$label<- paste ("Subsetted from--", y$label)

 return(ys)
}

