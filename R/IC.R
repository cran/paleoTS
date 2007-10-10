`IC` <-
function(w, logL=NULL, K=NULL, n=NULL, meth=c("AICc", "AIC", "BIC"))
# compute Information Criteria from log-likelihood, # parameters (K), and 
# sample size (n), if needed.
# can send just output of any opt.MODEL function (w) as well
{
 if(is.null(logL))
 	{	logL<- unname(w$val); K<- w$K; n<- w$n	}
 
 meth=match.arg(meth)	
 if (meth=="AIC")	ic<- -2*logL + 2*K
 if (meth=="AICc")	ic<- -2*logL + 2*K + (2*K*(K+1))/(n-K-1)
 if (meth=="BIC")	ic<- -2*logL + K*log(n)
 
 return(ic)
}

