# plotting functions #



#' Plot a paleoTS object
#'
#' @param x a \code{paleoTS} object
#' @param nse the number of standard errors represented by the error bars on the
#'   plot; defaults to 1
#' @param pool logical indicating if variances should be pooled across samples
#'   for the purposes of displaying error bars; defaults to \code{FALSE}
#' @param add logical, if \code{TRUE}, adds to existing plot
#' @param modelFit optional model fit from fitting functions
#' @param pch plotting symbol, defaults to 19
#' @param lwd line width, defaults to 1.5
#' @param ylim optional, y-limits of the plot
#' @param ... other arguments passed to plotting functions
#'
#' @return none.
#' @export
#' @importFrom graphics abline lines mtext plot points polygon segments text title
#'
#' @examples
#' x <- sim.GRW(ns = 30)
#' w <- fitSimple(x, model = "GRW", method = "Joint")
#' plot(x, modelFit = w)
plot.paleoTS<- function (x, nse = 1, pool = FALSE, add = FALSE, modelFit = NULL,
                         pch = 21, lwd = 1.5, ylim=NULL, ...)
{
  if (pool)
    x <- pool.var(x, ret.paleoTS = TRUE)
  se <- sqrt(x$vv/x$nn)
  lci <- x$mm - (nse * se)
  uci <- x$mm + (nse * se)
  xx <- x
  if (!is.null(x$start.age)) {
    if(x$timeDir=="decreasing") { x$tt <- x$start.age - x$tt; xl <- rev(range(x$tt))}
    else						{ x$tt<- x$tt + x$start.age; xl<- range(x$tt)}
  }
  else xl <- range(x$tt)
  if (!is.null(modelFit)) {
    mlab <- paste(modelFit$modelName, "expectation [95% prob. interval]")
    mc <- modelCurves(xx, w = modelFit)
    if (is.na(mc$ee[1]))
      modelFit <- NULL
  }
  if (is.null(modelFit))
    yl <- c(uci, lci)
  else yl <- c(uci, lci, mc$ll, mc$uu)
  if(is.null(ylim)) ylim<- range(yl)
  if (!add)
    plot(range(x$tt), ylim=ylim, typ = "n", pch = 19,
         xlab = "Time", ylab = "Trait Mean", xlim = xl, ...)
  if (!is.null(modelFit)) {
    if (!is.null(x$start.age)){
      if(x$timeDir == "decreasing") mc$tt <- x$start.age - mc$tt
      if(x$timeDir == "increasing") mc$tt <- x$start.age + mc$tt
    }
    polygon(c(mc$tt, rev(mc$tt)), c(mc$uu, rev(mc$ll)), col = grDevices::adjustcolor("grey", alpha.f = 0.4),
            border = "white")
    lines(mc$tt, mc$ee, col = "darkgrey", lwd = 2)
  }
  lines(x$tt, x$mm, lwd = lwd, ...)
  segments(x$tt, lci, x$tt, uci, lty = 1, lwd = lwd, ...)
  points(x$tt, x$mm, pch = pch, cex = 1.2, bg = "white", ...)
  mtext(x$label, cex = 0.7, col = "grey", font = 3)
  if (!is.null(modelFit))
    mtext(mlab, side = 4, cex = 0.8, col = "darkgrey", font = 2)
}


# internal function, not exported
modelCurves<- function(x, w, np=300)
  # returns list of model means, upper and lower 95% probability envelopes
{
  ee<- ii<- array(dim=np)  # set up arrays
  mn<- w$modelName
  mp<- w$par
  x0<- ifelse(w$method=="AD", x$mm[1], w$par["anc"])
  ttp<- seq(x$tt[1], x$tt[length(x$tt)], length.out=np)

  # list of models implemented
  okModels<- c("URW", "GRW", "Stasis", "OU", paste("Punc-", 1:100, sep=""))
  if (mn %in% okModels){
    if(mn=="URW"){ ee<- rep(x0, np); vv<- mp["vstep"]*ttp }
    if(mn=="GRW"){ ee<- x0+mp["mstep"]*ttp; vv<- mp["vstep"]*ttp }
    if(mn=="Stasis"){ ee<- rep(mp["theta"], np); vv<- rep(mp["omega"],np) }
    if(mn=="OU"){
      ee<- mp["theta"] * (1-exp(-mp["alpha"]*ttp)) + mp["anc"]*exp(-mp["alpha"]*ttp)
      vv<- (mp["vstep"]/(2*mp["alpha"])) * (1-exp(-2*mp["alpha"]*ttp)) }
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
