plot.paleoTS <-
function (x, nse = 1, pool = FALSE, add = FALSE, modelFit = NULL, 
    pch = 19, lwd = 1.5, ylim=NULL, ...) 
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
        if (!is.null(x$start.age)) 
            mc$tt <- x$start.age - mc$tt
        polygon(c(mc$tt, rev(mc$tt)), c(mc$uu, rev(mc$ll)), col = "wheat2", 
            border = "white")
        lines(mc$tt, mc$ee, col = "tan", lwd = 2)
    }
    lines(x$tt, x$mm, lwd = lwd, ...)
    segments(x$tt, lci, x$tt, uci, lty = 1, lwd = lwd, ...)
    points(x$tt, x$mm, pch = pch, cex = 1.2, ...)
    mtext(x$label, cex = 0.7, col = "grey", font = 3)
    if (!is.null(modelFit)) 
        mtext(mlab, side = 4, cex = 0.8, col = "tan", font = 2)
}
