as.paleoTSfit <-
function(logL, parameters, modelName, method, K, n, se)
{
  ic<- IC(logL=logL, K=K, n=n, method='AICc')
  y<- list(logL=logL, AICc=ic, parameters=parameters, modelName=modelName, method=method, K=K, n=n)
  class(y)<- "paleoTSfit"
  return(y)	
}
