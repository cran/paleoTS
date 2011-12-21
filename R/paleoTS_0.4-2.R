as.paleoTS<- function (mm, vv, nn, tt, MM=NULL, genpars=NULL, label=NULL, start.age=NULL, oldest="first", reset.time=FALSE)
# converts vectors of observed mean, variances, N's, ages (time steps),
# true means (for simulated series), and label to an object of class 'paleoTS'
{
  y<- list(mm=mm,vv=vv,nn=nn,tt=tt,MM=MM,genpars=genpars,label=label, start.age=start.age)
  
  if (oldest=="last")  # if samples are listed with youngest first, then reverse order
      {
      # reverse order of samples	
      oo<- length(y$mm):1
      y$mm <- y$mm[oo]
      y$vv <- y$vv[oo]
      y$nn <- y$nn[oo]
      y$tt <- y$tt[oo]	
    }
    
  if(reset.time){
  if (y$tt[1]!=0) {	# if starting sample not at t=0, make it so, and adjust start.age
	sa<- y$tt[1]
	if(!is.null(y$start.age) && sa!=y$start.age)	stop("Age of first sample does not match start.age")
	y$tt<- sa - y$tt
	y$start.age<- sa
	}
	}

  class(y)<- "paleoTS"
  return (y)
}


as.paleoTSfit<- function(logL, parameters, modelName, method, K, n, se)
{
  ic<- IC(logL=logL, K=K, n=n, method='AICc')
  y<- list(logL=logL, AICc=ic, parameters=parameters, modelName=modelName, method=method, K=K, n=n)
  class(y)<- "paleoTSfit"
  return(y)	
}


read.paleoTS<- function (file=NULL, oldest="first", reset.time=TRUE, ...)
# read in paleoTS data from a file
# samples should be listed from oldest (first) to youngest (last)
# if time is in ages BP (oldest with largest age), convert to forward moving time
#    by tt.new<- max(tt)-tt.old
{
  if (is.null(file))
     {
      ff<- file.choose()
      x<-read.table(ff, ...)
      lab1<- ff
     }

  else
     {
       x<-read.table(file=file, ...)
       lab1<-paste(getwd(), file)
     }
      
  # change to proper paleoTS object
  xr<- as.paleoTS(mm=x[,2], vv=x[,3], nn=x[,1], tt=x[,4], label=lab1, oldest=oldest, reset.time=TRUE)


  return (xr)
}

modelCurves<- function(x, w, np=500)
# returns list of model means, upper and lower 95% probability envelopes
{
  ee<- ii<- array(dim=np)  # set up arrays
  mn<- w$modelName
  mp<- w$par
  x0<- ifelse(w$method=="AD", x$mm[1], w$par["anc"])
  ttp<- seq(x$tt[1], x$tt[length(x$tt)], length.out=np)
  #ttp<- x$tt
  #ttp<- sort(c(seq(x$tt[1], x$tt[length(x$tt)], length.out=np), x$tt))
  #print(ttp)
  
  # list of models implemented
  okModels<- c("URW", "GRW", "Stasis", "OU", paste("Punc-", 1:100, sep=""))
  if (mn %in% okModels){  
	  if(mn=="URW"){ ee<- rep(x0, np); vv<- mp["vstep"]*ttp}  
	  if(mn=="GRW"){ ee<- x0+mp["mstep"]*ttp; vv<- mp["vstep"]*ttp}
	  if(mn=="Stasis"){ ee<- rep(mp["theta"], np); vv<- rep(mp["omega"],np) }
	  if(mn=="OU"){ if(!is.null(x$start.age)) tto<- x$start.age-ttp
	  				else	tto<- ttp
	  				ee<- mp["theta"] * (1-exp(-mp["alpha"]*tto)) + mp["anc"]*exp(-mp["alpha"]*tto)
	  				vv<- (mp["vstep"]/(2*mp["alpha"])) * (1-exp(-2*mp["alpha"]*tto)) }
	  if(grepl("Punc", mn)){
	  		# handle time shifts
	  		sp<- w$par[grep("shift", names(w$par))]  # shift samples
	  		st<- x$tt[sp]  # shift times
	  		ng<- length(sp)+1
	  		ggt<- cut(ttp, breaks=c(min(x$tt)-1, st, max(x$tt)+1), right=TRUE, labels=1:ng)
	
	
	  		# extract needed paramaters
	  		th<- w$par[grep("theta", names(w$par))]  # theta estimates
	  		om<- w$par[grep("omega", names(w$par))]  # omega estimates
	  		if(length(om)==1)	om<- rep(om, ng)	 # make into vector if needed
			
			# get ee and vv
			ee<- th[ggt]
			vv<- om[ggt]  		
		}  
	} else {ee<-NA; vv<- NA; warning(paste("modelFit argument not implemented for model", mn, "\n"))}
	
	
  if (!is.null(x$start.age))	tto<- x$start.age - ttp  else tto<- ttp 
  res<- list(tt=ttp, ee=ee, ll=ee-1.96*sqrt(vv), uu=ee+1.96*sqrt(vv))	
  return(res)	
}


plot.paleoTS<- function (x, nse=1, pool=FALSE, add=FALSE, modelFit=NULL, pch=19, lwd=1.5, ...)
# plots paleoTS object, with nse*se error bars
{
     
  if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
  se<- sqrt(x$vv/x$nn)
  lci<- x$mm-(nse*se)
  uci<- x$mm+(nse*se)
  xx<- x
  
  if (!is.null(x$start.age))     {x$tt<- x$start.age - x$tt; xl<- rev(range(x$tt))}
  else							xl<- range(x$tt)
  
  if (!is.null(modelFit))  
  {
   mlab<- paste(modelFit$modelName, "expectation [95% probability interval]")
   mc<- modelCurves(xx, w=modelFit)
   if(is.na(mc$ee[1]))	modelFit<- NULL   ## change back if modelFit not implemented for a model   	
  }
  		 
  if(is.null(modelFit)) yl<- c(uci, lci)
  else yl<- c(uci, lci, mc$ll, mc$uu)
 
  if(!add)	plot(range(x$tt), ylim=range(yl), typ="n", pch=19, xlab="Time", ylab="Trait Mean", xlim=xl, ...)

  if (!is.null(modelFit))  # add prob envelope, then expectation before data
   {
	if(!is.null(x$start.age))	mc$tt<- x$start.age - mc$tt
	polygon(c(mc$tt, rev(mc$tt)), c(mc$uu, rev(mc$ll)), col='wheat2', border="white")
	lines(mc$tt, mc$ee, col='tan', lwd=2)
   }
  lines (x$tt, x$mm, lwd=lwd, ...)
  segments (x$tt, lci, x$tt, uci, lty=1, lwd=lwd, ...)
  points (x$tt, x$mm, pch=pch, cex=1.2, ...)
  mtext(x$label, cex=0.7, col='grey', font=3)
  if(!is.null(modelFit)) mtext(mlab, side=4, cex=0.8, col='tan', font=2)

}




akaike.wts<- function(aa)
## aa is a vector of AIC or AICc values
{
 ma<- min(aa)	
 delt<- aa - ma		# delta values
 denom<- sum(exp(-delt/2))
 
 ww<- exp(-delt/2) / denom
 names(ww)<- names(aa)
 return(ww)	
}



IC<- function(logL, K, n=NULL, method=c("AICc", "AIC", "BIC"))
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



pool.var<- function(y,nn=NULL, ret.paleoTS=FALSE)
# calc pooled variance
# y is either a paleoTS object, or an array of variances
{
 # check if y is a paleoTS object
 if (class(y)=="paleoTS")
  {
    if (all(y$nn==1)) 	vp<- mean(y$vv)
    else    vp<- sum(y$vv*(y$nn-1))/sum(y$nn-1)
  }
 else
    vp<- sum(y*(nn-1))/sum(nn-1)
    
 if (ret.paleoTS)
 {
   yn<- y
   yn$vv<- rep(vp, length(y$mm))
   return(yn)		
 }
 else { return (vp)}
}

test.var.het<- function (y, method="Bartlett")
# test for variance heterogeneity among samples in a paleoTS object
{

 vp<- pool.var(y)
 NN<- sum(y$nn)
 k<- length(y$vv)
 top<- (NN-k)*log(vp)-(sum((y$nn-1)*log(y$vv)))
 bot<- 1+ (1/(3*(k-1)))*( (sum(1/(y$nn-1))) - (1/(NN-k)))
 TT<- top/bot
 p.val<- pchisq(TT, df=k-1, lower.tail=FALSE)

 w<-list(stat=TT, p.value=p.val, df=k-1)
 return (w)
}

ln.paleoTS <- function (y)
# returns paleoTS, with data approx ln-transformed
# mean(ln[y])= ln(mean[y]); var(ln[y])=(sd[y]/mean[y])^2
{
 logx<- y
 logx$mm<- log(y$mm)
 logx$vv<- (sqrt(y$vv)/y$mm)^2
 
 return (logx)
}

std.paleoTS <- function (y, zero="start")
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


sub.paleoTS <- function (y, ok=NULL, k=0.1)
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




compareModels<- function(..., silent=FALSE)
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


fit3models <- function (y, pool=TRUE, silent=FALSE, method=c("AD", "Joint"))
## fit stasis, GRW, URW models to a paleoTS
{
 method<- match.arg(method)
 
 if(method=="AD")
 {
  m1<- opt.GRW(y, pool=pool)
  m2<- opt.URW(y, pool=pool)
  m3<- opt.Stasis(y, pool=pool)	
 }
 	
 else if(method=="Joint")
 {
  m1<- opt.joint.GRW(y, pool=pool)
  m2<- opt.joint.URW(y, pool=pool)
  m3<- opt.joint.Stasis(y, pool=pool)	
 }
 
mc<- compareModels(m1,m2,m3, silent=silent)
invisible(mc)
}




logL.GRW<- function(p,y)
# function to return log-likelihood of step mean and variance (M,V)= p,
# given a paleoTS object
{
  # get parameter values: M is Mstep, V is Vstep
  M<- p[1]
  V<- p[2]
  
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  #S<- -0.5*log(2*pi*(V*dt+svAD)) - ((dy-(M*dt))^2)/(2*(V*dt + svAD))
  S<- dnorm(x=dy, mean=M*dt, sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}


logL.URW<- function(p,y)
# function to return log-likelihood of step mean and variance (M,V)= p,
# given a paleoTS object
{
  # get parameter values: V is Vstep
  V<- p[1]
    
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  #S<- -0.5*log(2*pi*(V*dt+svAD)) - (dy^2)/(2*(V*dt + svAD))
  S<- dnorm(x=dy, mean=rep(0,nd), sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}


logL.Stasis <- function(p, y)
## logL of stasis model
{
 # get parameter estimates
 M<-p[1]	# M is theta
 V<-p[2]	# V is omega
 dy<- diff(y$mm)
 nd<- length(dy)
 
 sv<- y$vv/y$nn
 svD<- sv[2:(nd+1)]  # only need sampling variance of descendant
 anc<- y$mm[1:nd]

 #S<- -0.5*log(2*pi*(V+svD)) - ((dy-(M-anc))^2)/(2*(V+svD))
 S<- dnorm(x=dy, mean=M-anc, sd=sqrt(V + svD), log=TRUE)
 return(sum(S))
}


logL.Mult<- function (p, y, model=c("GRW", "URW"))
# calculate logL over multiple sequences
#  here, y is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(y)
  for (i in 1:nseq)
    {
     if (model=="URW")
          Smult<- Smult + logL.URW(p,y[[i]])
     else if (model=="GRW")
          Smult<- Smult + logL.GRW (p,y[[i]])
    }
  return (Smult)
}


logL.SameMs <- function (p, y)
# computes logL over >1 sequence, of model in which all sequences have the 
# same directionality (Mstep), with different step variances
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m, v-1,..v-nseq}
{
  Sm<-0
  nseq<- length(y)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[1],p[i+1]), y[[i]])
   	
  return(Sm) 	
}


logL.SameVs <- function (p, y)
# computes logL over >1 sequence, of model in which all sequences have the 
# same step variance, with different steo means
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m-1,..m-nseq, vs}
{
  Sm<-0
  nseq<- length(y)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[i],p[nseq+1]), y[[i]])
   	
  return(Sm) 	
}





mle.GRW<- function(y)
# Gives analytical parameter estimates (GRW), assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(y$mm)-1 
 tt<- (y$tt[nn+1]-y$tt[1])/nn
 eps<- 2*pool.var(y)/round(mean(y$nn))  # sampling variance
 dy<- diff(y$mm)
 my<- mean(dy)
 
 mhat<- my/tt
 vhat<- (1/tt)*( (1/nn)*sum(dy^2) - my^2 - eps)
 
 w<- c(mhat, vhat)
 names(w)<- c("mstep", "vstep")
 return(w)
}


mle.URW<- function(y)
# Gives analytical parameter estimates (URW), assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(y$mm)-1 
 tt<- (y$tt[nn+1]-y$tt[1])/nn
 eps<- 2*pool.var(y)/round(mean(y$nn))  # sampling variance
 dy<- diff(y$mm)
 my<- mean(dy)
 
 vhat<- (1/tt)*( (1/nn)*sum(dy^2) - eps)
 
 w<- vhat
 names(w)<- "vstep"
 return(w)
}


mle.Stasis <- function (y)
# analytical solution to stasis model
{
  ns<- length(y$mm)
  vp<- pool.var(y)
  th<- mean(y$mm[2:ns])
  om<- var(y$mm[2:ns]) - vp/mean(y$nn)
  
  w<- c(th, om)
  names(w)<- c("theta", "omega")
  return(w)
}



opt.GRW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.GRW(y)
  if (p0[2] <= 0)	p0[2]<- 1e-7
  names(p0)<- c("mstep", "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.GRW failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}


opt.URW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.URW(y)
  if (p0 <= 0)	p0<- 1e-7
  names(p0)<- "vstep"
  if (is.null(cl$ndeps))		cl$ndeps<- p0/1e4
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.URW failed ", immediate.=TRUE)
		  w$par<- NA
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='AD', K=1, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}



opt.Stasis<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.Stasis(y)
  if (p0[2] <= 0 || is.na(p0[2]))	p0[2]<- 1e-7
  names(p0)<- c("theta", "omega")
  if (is.null(cl$ndeps))		cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.Stasis failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}



opt.RW.Mult<- function (yl, cl=list(fnscale=-1), model=c("GRW", "URW"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates single model across multiple sequences
# pool=TRUE will pool variances _within_ sequences
{
  if (class(yl)=="paleoTS")
     stop("opt.RW.mult is onlt for multiple paleoTS sequences\n")
  nseq<- length(yl)
     
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
     	ll<- c(NA,0)
     	p0<- mle.GRW(yl[[1]])	
     	K<- 2	 	}
  else if (model=="URW")
     { ll<- 0
       p0<- mle.URW(yl[[1]])
       K<- 1	  	}
  
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, y=yl, model=model), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, y=yl, model=model), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, y=yl, model=model), silent=TRUE)
    else 
      w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, y=yl, model=model), silent=TRUE)
	#if (class(w)=="try-error")  # if still doesn't work
	#  	{   warning("opt.RW.Mult failed ", immediate.=TRUE)
	#	    w$par<- NA
	#	    w$value<- NA }
  }
  
  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(model,'.Mult', sep=''), method='AD', K=K, n=n, se=w$se)
  
  return (wc)

  
}


opt.RW.SameMs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Ms model across multiple sequences
{
  if (class(yl)=="paleoTS")
  	stop("Function opt.SameMs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(yl[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(median(p0m), p0v)  # shared Ms, followed by separate Vs for each sequence
  names(p0)<- c("mstep", paste("vstep", 1:nseq, sep=""))
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)

  # optimize logL
  ll<- c(NA, rep(0,nseq))
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
    else 
      w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)
	#if (class(w)=="try-error")  # if still doesn't work
	#  	{   warning("opt.RW.SameMs failed ", immediate.=TRUE)
	#	    w$par<- NA
	#	    w$value<- NA }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameMs.Mult", method='AD', K=nseq+1, n=n, se=w$se)
  
  return (wc)
}


opt.RW.SameVs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Vs model across multiple sequences
{
  if (class(yl)=="paleoTS")
  	stop("Function opt.SameVs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(yl[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(p0m, median(p0v))  # separate Ms, followed by shared Vs for each sequence
  names(p0)<- c(paste("mstep", 1:nseq, sep=""), "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)  

  # optimize logL
  ll<- c(rep(NA,nseq), 0)
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)

  # if optim fails, set ndeps based on p0
  if (class(w)=="try-error")
  {
    cl$ndeps<- p0/10000
    if (meth=="L-BFGS-B")
      w<- try(optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, y=yl), silent=TRUE)
    else
      w<- try(optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, y=yl), silent=TRUE)
	#if (class(w)=="try-error")  # if still doesn't work
	#  	{   warning("opt.RW.SameVs failed ", immediate.=TRUE)
	#	    w$par<- NA
	#	    w$value<- NA }
  }

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameVs.Mult", method='AD', K=nseq+1, n=n, se=w$se)

  return (wc)
}



sim.GRW <- function (ns=20, ms=0, vs=0.1, vp=1, nn=rep(20,ns), tt=0:(ns-1))
# simulates GRW; ns= number of samples, ms=mean and vs=variance of the step distribution,
#  vp=population variance, tt=ages of the samples
{
 MM<- array(dim=ns)
 mm<- array(dim=ns)
 vv<- array(dim=ns)
 dt<- diff(tt)
 
 inc<- rnorm(ns-1, ms*dt, sqrt(vs*dt))	# evolutionary increments
  
 MM<- cumsum(c(0,inc))	# true means
 mm<- MM + rnorm(ns, 0, sqrt(vp/nn))	# true means plus sampling error
 vv<- rep(vp, ns)
 
 gp<- c(ms, vs)
 names(gp)<- c("mstep", "vstep")
 
 res<- as.paleoTS(mm=mm, vv=vv, nn=nn, tt=tt, MM=MM, genpars=gp, label="Created by sim.GRW()") 
 return(res)
}


sim.Stasis <- function(ns=20, theta=0, omega=0, vp=1, nn=rep(20,ns), tt=0:(ns-1))
# simulate stasis
{
 xmu<- rnorm(ns, mean=theta, sd=sqrt(omega))
 xobs<- xmu + rnorm(ns, 0, sqrt(vp/nn))
 gp<- c(theta, omega)
 names(gp)<- c("theta", "omega")
   	
 x <- as.paleoTS(mm=xobs,vv=rep(vp,ns),nn=nn,tt=tt,MM=xmu,genpars=gp,label="Created by sim.stasis()") 
 return(x)	 	 	
}

lynchD<- function (y, gen.per.t=1e6, pool=TRUE, method=c('AD', 'Joint'), ...)
# compute Lynch's rate metric, Delta
# gen.per.t is the number of generations per unit tt
# ... are further arguments to opt.URW()
{
  method<- match.arg(method)
  vp<- pool.var(y)
  if (method == 'AD')
    {
      wu<- opt.URW(y, pool=pool)
      vs<- unname(wu$par)	
    }
  if (method == 'Joint')
    {
      wu<- opt.joint.URW(y, pool=pool)
      vs<- unname(wu$par[2])	
    }
    
  D<- 0.5*(vs/vp)/gen.per.t
  
  drift.min<- 5e-5
  drift.max<- 5e-3
  
  if (D < drift.min)	res<- "Slower than drift range"
  else if (D > drift.max)	res<- "Faster than drift range"
  else					res<- "Within range of drift"
  
  w<- list(D=D, pooled.var=vp, gen.per.t=gen.per.t, vstep=vs, drift.range=c(drift.min, drift.max), result=res) 
  return(w)	
}

## functions added from paleoTSalt


logL.joint.GRW<- function (p, x)
# returns logL of GRW model for paleoTS object x
# p is vector of parameters: alpha, ms, vs
{
 # prepare calculations
 anc<- p[1]
 ms<- p[2]
 vs<- p[3]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- vs*outer(x$tt, x$tt, FUN=pmin)
 diag(VV)<- diag(VV) + x$vv/x$nn 
 
 # compute logL based on multivariate normal
 M<- rep(anc, n) + ms*x$tt
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
 
 return(S)		
}

logL.joint.URW<- function(p,x)
# returns logL of URW model for paleoTS object x
# p is vector of parameters: alpha, vs
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- vs*outer(x$tt, x$tt, FUN=pmin)
 diag(VV)<- diag(VV) + x$vv/x$nn 	# add sampling variance
 
 # compute logL based on multivariate normal
 M<- rep(anc, n)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)

 return(S)
}

logL.joint.Stasis<- function (p,x)
# returns logL of Stasis model for paleoTS object x
# p is vector of parameters: theta, omega
{
 # prepare calculations
 theta<- p[1]
 omega<- p[2]
 n<- length(x$mm)
 
 # compute covariance matrix
 VV<- diag(omega + x$vv/x$nn)  # omega + sampling variance
 
 # compute logL based on multivariate normal
 M<- rep(theta, n)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
 return(S)	 	
}


opt.joint.GRW<- function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize GRW model using alternate formulation
{
 ## check if pooled, make start at tt=0
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt)
  
 ## get initial estimates
 p0<- array(dim=3)
 p0[1]<- x$mm[1]	
 p0[2:3]<- mle.GRW(x)
 if (p0[3]<=0)	p0[3]<- 1e-7
 names(p0)<- c("anc", "mstep", "vstep")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
 if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, lower=c(NA,NA,0), hessian=hess, x=x)
 else 					w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, hessian=hess, x=x)

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='Joint', K=3, n=length(x$mm), se=w$se)
  
 return (wc)
}

opt.joint.URW<- function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize URW model using alternate formulation
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt)
 
 ## get initial estimates
 p0<- array(dim=2)
 p0[1]<- x$mm[1]
 p0[2]<- mle.URW(x)
 names(p0)<- c("anc","vstep")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
 if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, lower=c(NA,0), hessian=hess, x=x)
 else					w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, hessian=hess, x=x)

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='Joint', K=2, n=length(x$mm), se=w$se)
  
 return (wc)
}

opt.joint.Stasis<- function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize Stasis model using alternate formulation
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt)
  
 ## get initial estimates
 p0<- mle.Stasis(x)
 if(p0[2]<=0 || is.na(p0[2]))	p0[2]<- 1e-7
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-9
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, lower=c(NA,0), hessian=hess, x=x)
 else				  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, hessian=hess, x=x)
 
 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='Joint', K=2, n=length(x$mm), se=w$se)
  
 return (wc)
	
}


 # functions to compute mean and variance of OU process
 ou.M<- function(anc, theta, aa, tt) theta*(1 - exp(-aa*tt)) + anc*exp(-aa*tt)
 ou.V<- function(vs, aa, tt)        (vs/(2*aa))*(1 - exp(-2*aa*tt))


sim.OU<- function (ns=20, anc=0, theta=10, alpha=0.3, vs=0.1, vp=1, nn=rep(20, ns), tt=0:(ns-1))
## generate a paleoTS sequence according to an OU model
{
    
    MM <- array(dim = ns)
    mm <- array(dim = ns)
    vv <- array(dim = ns)
    dt <- diff(tt)
    MM[1] <- anc
    x <- rnorm(nn[1], mean = MM[1], sd = sqrt(vp))
    mm[1] <- mean(x)
    vv[1] <- var(x)
    for (i in 2:ns) {
        ex<- ou.M(MM[i-1], theta, alpha, dt[i-1])
        vx<- ou.V(vs, alpha, dt[i-1])
        MM[i]<- rnorm(1, ex, sqrt(vx))    
        x <- rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
        mm[i] <- mean(x)
        vv[i] <- var(x)
    }
    
    gp <- c(anc, theta, alpha, vs)
    names(gp) <- c("anc", "theta", "alpha", "vs")
    res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
        genpars = gp, label = "Created by sim.OU()")

    return(res)
}

logL.joint.OU<- function(p,x)
# returns logL of OU model for paleoTS object x
# 
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 theta<- p[3]
 aa<- p[4]
 n<- length(x$mm)
 
 # compute covariance matrix
 ff<- function (a,b) abs(a-b)
 VV<- outer(x$tt, x$tt, FUN=ff)
 VV<- exp(-aa*VV)
 VVd<- ou.V(vs,aa,x$tt)
 VV2<- outer(VVd,VVd,pmin)
 VV<- VV*VV2
 diag(VV)<- VVd + x$vv/x$nn 	# add sampling variance
 #detV<- det(VV)
 #invV<- solve(VV)
 
 # compute logL based on multivariate normal
 M<- ou.M(anc, theta, aa, x$tt)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 S<- dmvnorm(t(x$mm), mean=M, sigma=VV, log=TRUE)
 return(S)		
}


opt.joint.OU<- function (x, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize OU model using tree methods
{
 ## check if pooled
 if (pool)	x<- pool.var(x, ret.paleoTS=TRUE)
 x$tt<- x$tt - min(x$tt) 
 
 ## get initial estimates
 w0<- mle.GRW(x)
 halft<- (x$tt[length(x$tt)]-x$tt[1])/4			# set half life to 1/4 of length of sequence
 p0<- c(x$mm[1], w0[2]/10, x$mm[length(x$mm)], log(2)/halft)
 names(p0)<- c("anc","vstep","theta","alpha")
 #print(p0)
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, lower=c(NA,1e-10,NA,1e-8), hessian=hess, x=x)
 else 				  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, hessian=hess, x=x) 

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='OU', method='Joint', K=4, n=length(x$mm), se=w$se)
  
 return (wc)
}


##### start of punctuate functions

cat.paleoTS<- function (y)
# concatenates multiple paleoTS objects, with y a list of paleoTS objects
{
 x<- y[[1]]
 for (i in 2:length(y))
   {
   	x$mm<- append(x$mm, y[[i]]$mm)
   	x$vv<- append(x$vv, y[[i]]$vv)
   	x$tt<- append(x$tt, y[[i]]$tt)
   	x$MM<- append(x$MM, y[[i]]$MM)
   	x$nn<- append(x$nn, y[[i]]$nn)
   }
  
  return (x)   	
}


sim.punc<- function (ns=c(10,10), theta=c(0,1), omega=rep(0,length(theta)), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate punctuated sequence; theta and omega are vectors of paramters
# ns is vector of ns in each sub-sequence
{
  nr<- length(theta)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1		
   	 }
   	 
   	 xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta[i], omega=omega[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.punc()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(theta, omega, shft)
  names(y$genpars)<- c(paste("theta",1:nr,sep=""), paste("omega",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
  return (y)  	
}

split4punc<- function (y, gg, overlap=TRUE)
# divides a paleoTS object (y) into a several paleoTS objects, according to vector 'gg'
# gg is a vectors of 1,2,3.. indicating groupings
# overlap=TRUE means that the adjacent samples are included 
{
  yl<- list()
  ng<- max(gg)
  for (i in 1:ng)
   {
   	 ok<- gg==i
   	 # it's counterintuitive, but the next line is correct given the interface with other functions
   	 if (i<ng & overlap==TRUE)   ok[min(which(gg==i+1))]<- TRUE  # right!
   	 #if(i>1 & overlap==TRUE)	ok[max(which(gg==i-1))]<- TRUE   # this is not right!
   	
   	 yl[[i]]<- as.paleoTS(y$mm[ok],y$vv[ok],y$nn[ok],y$tt[ok],y$MM[ok])
   }
  return (yl)	
}


logL.punc<- function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  xl<- split4punc(y,gg)
  print(xl[[1]]$mm)
  S<-0
  for (i in 1:ng)
    S<- S+ logL.Stasis(p=c(th[i], om[i]), xl[[i]])
  
  return (S)	
}

logL.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S + logL.Stasis(p=c(th[i], om), xl[[i]])
  
  return (S)	  	
}

opt.punc<- function(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare) 
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
 
  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- 2*ng; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 3*ng-1; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn
  
  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }
  
  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)

  return(wc)
}

logL.joint.punc<- function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  M<- th[gg]  # vector of MVN means
  VV<- diag(om[gg] + y$vv/y$nn)  # vcv matrix
  S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  
  return (S)	
}

logL.joint.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  M<- th[gg]  # vector of MVN means
  omv<- rep(om, max(gg))
  VV<- diag(omv[gg] + y$vv/y$nn)  # vcv matrix
  S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  
  return (S)	  	
}


opt.joint.punc<- function(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare) 
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
 
  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- 2*ng; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 3*ng-1; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn
  
  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.joint.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.joint.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }
  
  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)
}


shifts<- function (ns, ng, minb=5)
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

shift2gg<- function (ss, ns)
# ss is vector of shift points, ns is # samples
{
  z<- c(0,ss,ns+1)
  cc<- cut(1:ns, breaks=z, right=FALSE)
  gg<- as.numeric(cc)
  return(gg)	
}

fitGpunc<- function(y, ng=2, minb=5, pool=TRUE, oshare=TRUE, method=c('AD', 'Joint'), silent=FALSE, hess=FALSE, ...)
## optimize punctuation models (with some min n per section)
{
 method<- match.arg(method) 
 
 if(ng==1)  # if only one grouping, same as stasis model
 {
   warning('Fitting stasis model (because ng=1)')
   if (method=='AD')	ww<- opt.Stasis(y, pool=pool, hess=hess, ...)
   else if (method=='Joint')	ww<- opt.joint.Stasis(y, pool=pool, hess=hess, ...)
   return(ww)	
 }
 
 ns<- length(y$mm)
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("Total # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    # the different gg for AD and J is required for the interpretation of the "shift" parameters to be the same across parameterizations
    ggA<- shift2gg(GG[,i], ns)
    ggJ<- shift2gg(GG[,i]+1, ns)
    if(method=='AD')    		 w<- opt.punc(y, ggA, oshare=oshare, pool=pool, hess=hess, ...)
    else if (method=='Joint')    w<- opt.joint.punc(y, ggJ, oshare=oshare, pool=pool, hess=hess, ...)
    logl[i]<- w$logL
    wl[[i]]<- w
  }
 cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}

sim.GRW.shift <- function (ns=c(10,10), ms=c(0,1), vs=c(0.5,0.5), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate sequence of GRW with different parameter values in different segments
# ns is vector of ns in each sub-sequence, similar for other parameters
{
  nr<- length(ms)
  cns<- cumsum(ns)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- cns[i-1]+1
   	  end.i<- cns[i-1]+ns[i]	
   	 }
   	 
   	 xl[[i]]<- sim.GRW(ns=ns[i], ms=ms[i], vs=vs[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   	 if (i>1)
   	 	xl[[i]]$mm<- xl[[i]]$mm + xl[[i-1]]$mm[length(xl[[i-1]]$mm)]
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.GRW.shift()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(ms, vs, shft)
  names(y$genpars)<- c(paste("mstep",1:nr,sep=""), paste("vstep",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))  
  return (y)  	
}

opt.GRW.shift<- function(y, ng=2, minb=5, model=1, pool=TRUE, silent=FALSE)
## optimize for shifted GRW dynamics (with some min n per section)
## models:	1  grw (same Vs, diff Ms)
#			2  grw (same Ms, diff Vs)
#			3  urw (diff Vs)
#			4  grw (diff Ms, diff Vs)
{
 ns<- length(y$mm)
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    gg<- shift2gg(GG[,i], ns)
    yl<- split4punc(y, gg)
    #print(yl[[1]])
    if (model==1)	{ w<- opt.RW.SameVs(yl, pool=pool); w$modelName<- paste('GRWsameVs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model==2)	{ w<- opt.RW.SameMs(yl, pool=pool); w$modelName<- paste('GRWsameMs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model>=3)
     {
      wli<- list()
      totS<-0
      totpar<- numeric()
      for (j in 1:ng) 
       {
      	if (model==3){	
      		wli[[j]]<- opt.URW(yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste("vstep", j, sep="")   }
      	if (model==4){	
      		wli[[j]]<- opt.GRW (yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste(c("mstep","vstep"), j, sep="")   }
      	totS<- totS + wli[[j]]$logL
      	totpar<- c(totpar, wli[[j]]$parameters)
      	kk<- length(totpar)+ng-1
       }
    
      ifelse(model==3, mn<- 'URW-shift', mn<- 'GRW-shift')
      w<- as.paleoTSfit(logL=totS, parameters=totpar, modelName=paste(mn, ng-1, sep='-'), method='AD', K=kk, n=ns-1, se=NULL)	
     }
         
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
 
}

sim.sgs <- function (ns=c(20,20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate stasis-grw-stasis sequence, take theta2 to be final value after grw part
{
  xl<- list()
  for (i in 1:3)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1	
   	 }
   	 
   	 if (i==2)
   	 	xl[[i]]<- sim.GRW(ns[2], ms, vs, nn=nn[start.i:end.i], tt=tt[start.i:end.i], vp=vp)
   	 
   	 else
   	 	xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta, omega=omega, vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }

  ## add offsets
  xl[[2]]$mm<- xl[[2]]$mm + xl[[1]]$MM[ns[1]]
  xl[[3]]$mm<- xl[[3]]$mm + xl[[2]]$MM[ns[2]]
	
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.sgs()"
  y$genpars <- c(theta, omega, ms, vs)
  names(y$genpars)<- c("theta","omega", "ms","vs") 
  return(y)
}


logL.sgs<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1] 
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5:6] }
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4:5]	}
  
  l1<- logL.Stasis(p=c(th[1], om[1]), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])  
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om[2]), yl[[3]])
  
  logl<- l1+l2+l3
  return(logl)  	
}

logL.sgs.omega<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model, omega shared over stasis segments
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1]
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5]	}
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4]	}
  
  l1<- logL.Stasis(p=c(th[1], om), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om), yl[[3]])
  
  logl<- l1+l2+l3
  return(logl)  		
}


opt.sgs<- function(y,gg,cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare=TRUE, model="GRW")
# do optimization for sgs model
{
 yl<- split4punc(y, gg)
 hat<- mle.GRW(yl[[2]])
 if (hat[1] == 0)	hat[1]<- 1e-3
 if (hat[2] < 1e-4)	hat[2]<- 1e-4
 if (model=="URW")	p0<- c(hat[2], mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 else 				p0<- c(hat, mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 
 if (oshare)	
 	{ p0<- append(p0, mean(var(yl[[1]]$mm), var(yl[[3]]$mm)))
 	  if (model=="GRW")		{ K<-7; pn<- c("ms","vs","theta1","theta2","omega"); lw<- c(NA,0,NA,NA,0) }
 	  else					{ K<-6; pn<- c("vs","theta1","theta2","omega"); lw<- c(0,NA,NA,0) }	}  
 else
 	{ p0<- append(p0, c(var(yl[[1]]$mm), var(yl[[3]]$mm)) )
 	  if (model=="GRW")	{ K<- 8; pn<- c("ms","vs","theta1","theta2","omega1","omega2"); lw<- c(NA,0,NA,NA,0,0) }
 	  else 				{ K<- 7; pn<- c("vs","theta1","theta2","omega1","omega2"); lw<- c(0,NA,NA,0,0) }
 	}
 
 cl$ndeps <- p0/100
 names(p0)<- pn

 if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.sgs.omega, gg=gg, method=meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn=logL.sgs.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.sgs, gg=gg, method= meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn = logL.sgs, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }

  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else 			w$se<- NULL

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('SGS', model, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)  
  return(wc) 
}

fit.sgs<- function(y, minb=5, oshare=TRUE, pool=TRUE, silent=FALSE, hess=FALSE, meth="L-BFGS-B", model="GRW")
## optimize for stasis-GRW-stasis dynamics (with some min n per section)
{
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)  # pool variances
 ns<- length(y$mm)
 ng<-3
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n\n", "i\tshifts\tlogL\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc) 
  {
    gg<- shift2gg(GG[,i], ns)
	w<- opt.sgs(y, gg, oshare=oshare, hess=hess, meth=meth, model=model)              
    if (!silent) 	cat (i, "\t", GG[,i], "\t", round(w$logL, 3), "\n")
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 if (!silent) cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}


sim.covTrack<- function (ns=20, b=1, evar=0.1, z, nn=rep(20, times=ns), tt=0:(ns-1), vp=1)
# simulates tracking optimum, with covariate z
{
  # check on length of covariate
  if(length(z)==ns)	
  	{ 
  	 z<- diff(z)
	 warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
  
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- rnorm(ns-1, b*z, sqrt(evar))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))  # add sampling error
  vv <- rep(vp, ns)
  gp <- c(b, evar)
  names(gp) <- c("slope.b", "evar")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
        genpars = gp, label = "Created by sim.covTrack()")
  return(res)	
}




logL.covTrack<- function(p, y, z)
# z is covariate; pars = b (slope), evar (variance)
# IMPT: z is of length ns-1; one for each AD transition IN ORDER
{
  b<- p[1]
  evar<- p[2]
  dy <- diff(y$mm)
  dt <- diff(y$tt)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svA <- sv[1:nd]
  svD <- sv[2:(nd + 1)]
  svAD <- svA + svD
  #S <- -0.5 * log(2 * pi * (evar + svAD)) - ((dy - (b*z))^2)/(2 * (evar + svAD))
  S<- dnorm(x=dy, mean=b*z, sd=sqrt(evar+svAD), log=TRUE)
  return(sum(S))	
}


opt.covTrack<- function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE) 
{
    # check if z is of proper length; first difference if necessary
    ns<- length(y$mm)
    if(length(z)==length(y$mm))	
  	 { 
  	  z<- diff(z)
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(z) != ns-1)  stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n" )

    
    # get initial estimates by regression
    reg<- lm(diff(y$mm) ~ z-1)
    p0<- c(coef(reg), var(resid(reg)))
    names(p0) <- c("b", "evar")
    
    # pool variances if needed and do optimization
    if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
    if (is.null(cl$ndeps)) 
        cl$ndeps <- abs(p0/10000)
    if (meth == "L-BFGS-B") 
        w <- optim(p0, fn=logL.covTrack, method = meth, lower = c(NA, 0), control = cl, hessian = hess, y=y, z=z)
    else  w<- optim(p0, fn=logL.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)
   

    if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
    else w$se <- NULL
    
	wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='AD', K=2, n=length(y$mm)-1, se=w$se)
}

logL.Mult.covTrack<- function (p, yl, zl)
# y is a _list_ of paleoTS objects, and z is a _list_ of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
    {	Smult<- Smult + logL.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)	
}


opt.covTrack.Mult<- function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
# y and z are lists of paleoTS, and covariates, respectively
{
 if (class(yl) == "paleoTS") 
      { stop("opt.track.Mult is only for multiple paleoTS sequences\n") }
 nseq <- length(yl)
 if (pool) {
        for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE) }
        
 # check lengths of z
 for (i in 1:nseq)
 {
    ns<- length(yl[[i]]$mm)
    if(length(zl[[i]])==length(yl[[i]]$mm))	
  	 { 
  	  zl[[i]]<- diff(zl[[i]])
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(zl[[i]]) != ns-1)  stop("Covariate length [", length(zl[[i]]), "] does not match the sequence length [", ns, "]\n" )
 	
 }

        
 p0<- c(0,var(yl[[1]]$mm))  # lousy method for initial guess!
 names(p0)<- c("b", "evar")
 w<- optim(p0, fn=logL.Mult.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,0), hessian=hess, yl=yl, zl=zl)
 
 if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
 else		w$se<- NULL
 
 ff<- function(x) length(x$mm)-1
 n<- sum(sapply(yl, ff))
 
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='trackCovariate.Mult', method='AD', K=2, n=n, se=w$se)
 return(wc)	
}


LRI<- function(x, gen.per.t=1e6, draw=TRUE)
{
  n<- length(x$mm)
  x$tt<- x$tt*gen.per.t  # convert to generational timescale (e.g., 1e6 when time is measured in Myr and generation time = 1 yr)
  sp<- sqrt(pool.var(x))
  dy<- outer(x$mm, x$mm, FUN='-')  # outer() much faster than looping
  dt<- outer(x$tt, x$tt, FUN='-')
  dy<- abs(dy)
  dt<- abs(dt)
   
  dy<- dy/sp
  ut<- upper.tri(dy)
  logI<- log10(dt[ut])
  logR<- log10(dy[ut]) - logI
   
  # subset only non-zero rates
  ok<- is.finite(logR)
  if (sum(ok) != length(logI))	warning("Some zero rates were ignored.")
  
  # function to fit line in LRI plot robustly, using min abs deviations as per Gingerich 1993
  opt.lad<- function (x, y)
	{
	  ok<- is.finite(x) & is.finite(y)
	  xok<- x[ok]
	  yok<- y[ok]
	  w.ls<- lm(yok ~ xok)  # use LS coef as starting point in optimization
	
	  # function minimized for least abs deviation
	  fad<- function(p, x, y)  sum( abs(y - p[2]*x - p[1]) )
	  w.lad<- optim(w.ls$coef, fn=fad, x=xok, y=yok)
	  return(w.lad$par)	
	}
    
  # call robust regression function
  lf<- opt.lad(logI[ok], logR[ok])
  lf[3]<- 10^lf[1]   # generational rate
  names(lf)<- c("Intercept", "slope", "GenerationalRate")
  
  # do LRI plot, if desired
  if (draw)
   {
   	xl<- c(0, max(logI[ok]))
   	yl<- c(min(logR[ok]), max(lf[1], max(logR[ok])))
   	plot(logI[ok], logR[ok], xlim=xl, ylim=yl, xlab='log10 Interval [generations]', ylab='log10 Rate', cex=0.6)
   	abline(lf[1], lf[2], col='black', lwd=2)
   	title("LRI plot")
   	mtext(paste('data label: ', x$lab), cex=0.6, font=3, col='darkgrey')
   	restext<- paste('Slope = ', round(lf[2],3), '\n', 'Intercept = ', round(lf[1],3), '\n', 'Generational Rate = ', round(lf[3],5), '\n', sep='')
   	text(0, min(logR[ok]), restext, adj=c(0,0), cex=0.7, font=2)
   }
  
  #w<- list(b0=b0, b1=b1, h0=h0, logR=logR, logI=logI, dy=dy)  
  return(lf)
}

