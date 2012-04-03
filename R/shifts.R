shifts <-
function (ns, ng, minb=5)
## returns matrix of all allowable shift points for 
## dividing ns samples into ng groups of at least minb samples each
{
  aa<- combn(ns, ng-1)
  #if (ng==2)	ok<- aa>=minb & aa<=ns-minb+1
  if (ng==2)	ok<- aa>minb & aa<=ns-minb+1
  if (ng>2)
  {
    daa<- apply(aa,2,diff)
    if (ng>3)	mdaa<- apply(daa,2,min)
    else 		mdaa<- daa
    #ok1<- mdaa >= minb			# only long enough combinations
    ok1<- mdaa >= minb
    ok2<- aa[1,] > minb 		# must start far enough from 1 for minb
    #ok2<- aa[1,] >= minb
    #ok3<- aa[ng-1, ] <= ns-minb+1		# must end short of end
    ok3<- aa[ng-1,]  <= ns-minb+1
    ok<- ok1&ok2&ok3
   }
  ret<- aa[,ok]
  if (ng==2)	ret<- matrix(ret, nrow=1)	# convert to matrix from vector is ng=2
  if (is.null(dim(ret)))	ret<- matrix(ret, ncol=1) # handle if only 1 possible grouping and ng>2
  if (sum(ok)>0)	return(ret)
  else 				return(NULL)
}
