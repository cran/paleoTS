IC <-
function(logL, K, n=NULL, method=c("AICc", "AIC", "BIC"))
# compute Information Criteria from log-likelihood, # parameters (K), and 
# sample size (n), if needed.
{
 method<-match.arg(method)	
 if ((method=="AICc" || method=="BIC") && is.null(n))  stop('AICc requires n.')
 if (method=="AIC")	ic<- -2*logL + 2*K
 if (method=="AICc")	ic<- -2*logL + 2*K + (2*K*(K+1))/(n-K-1)
 if (method=="BIC")	ic<- -2*logL + K*log(n)
 
 return(ic)
}
