compareModels <-
function(..., silent=FALSE)
{
  modelList<- list(...)
  
  # make sure all are paleoTSfit objects, and all use same method (AD or Joint)
  classv<- sapply(modelList, FUN=class)
  methv<- sapply(modelList, FUN=function(x){x$method})
  nv<- sapply(modelList, FUN=function(x){x$n})
  nm<- length(modelList)
  
  if(!all(classv=='paleoTSfit'))  	stop("All objects must be of class 'paleoTSfit'")
  if(!all(methv==methv[1]))			stop(paste("All model fits must use the same method (AD or Joint)", sep='\n'))
  else method<- methv[1]
  if(!all(nv==nv[1]))				stop("Objects have differing n.")
  else nn<- nv[1]
 
  
  # construct data frame and parameter list
  logL<- sapply(modelList, FUN=function(x){x$logL})
  K<- sapply(modelList, FUN=function(x){x$K})
  AICc<- sapply(modelList, FUN=function(x){x$AICc})
  Akaike.wt<- round(akaike.wts(AICc),3)
  df<- data.frame(logL, K, AICc, Akaike.wt)
  row.names(df)<- sapply(modelList, FUN=function(x){x$modelName})
  
  pl<- lapply(modelList, FUN=function(x){x$parameters})
  names(pl)<- row.names(df)
 
  # print information
    if(!silent)
  	{
  	  	cat ('\nComparing ', nm, ' models [n = ', nn, ',', ' method = ', methv[1], ']\n\n', sep='')
  	  	print (df)
  	}
 
 if(silent)		return(list(modelFits=df, parameters=pl))
 else 			invisible(df)
}
