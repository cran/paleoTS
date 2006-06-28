`fit3models` <-
function (y, pool=TRUE, silent=FALSE, wts="AICc")
## fit stasis, GRW, URW models to a paleoTS
{
 # fit RW, RWu, stasis models
 m.grw<- opt.RW(y, pool=pool)
 m.urw<- opt.RWu(y, pool=pool)
 m.st <- opt.stasis(y, pool=pool)
 
 # collect information to return
 aic<- c(m.grw$AIC, m.urw$AIC, m.st$AIC)
 aicc<- c(m.grw$AICc, m.urw$AICc, m.st$AICc)
 logl<- c(m.grw$value, m.urw$value, m.st$value)
 hats<- c(m.grw$par, m.urw$par, m.st$par) 
 if (wts=="AICc")	ak.wts<- akaike.wts(aicc)
 else				ak.wts<- akaike.wts(aic)
 
 w<- list(aic=aic, aicc=aicc, logl=logl, hats=hats, ak.wts=ak.wts)
 if (!silent)
  {
  cat ("Results Summary:\n\n")
  rt<- cbind(logl, aic, aicc, ak.wts)
  row.names(rt)<- c("GRW", "URW", "Stasis")
  print(rt)
  hh<- hats
  names(hh)<- c("Ms", "Vs", "Vsu", "Theta", "Omega")
  cat ("\n\nParameter estimates: \n")	
  print (hh)
  }
 
 else
 	return(w)
}

