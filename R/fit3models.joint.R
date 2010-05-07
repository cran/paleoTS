fit3models.joint <-
function (y, pool = TRUE, silent = FALSE, wts = "AICc") 
{
    mn<- c("GRW", "URW", "Stasis")
    m.grw <- opt.joint.GRW(y, pool = pool)
    m.urw <- opt.joint.URW(y, pool = pool)
    m.st <- opt.joint.Stasis(y, pool = pool)
    aic <- c(m.grw$AIC, m.urw$AIC, m.st$AIC)
    aicc <- c(m.grw$AICc, m.urw$AICc, m.st$AICc)
    logl <- c(m.grw$value, m.urw$value, m.st$value)
    hats <- c(m.grw$par, m.urw$par, m.st$par)
    names(hats)<- c('anc.GRW', 'mstep.GRW', 'vstep.GRW', 'anc.URW', 'vstep.URW', 'theta.Stasis', 'omega.Stasis')
    if (wts == "AICc") 
        ak.wts <- akaike.wts(aicc)
    else ak.wts <- akaike.wts(aic)
    names(aic)<- names(aicc)<- names(logl)<- names(ak.wts)<- mn

    w <- list(aic = aic, aicc = aicc, logl = logl, hats = hats, 
        ak.wts = ak.wts)
    if (!silent) {
        cat("Results Summary:\n\n")
        rt <- cbind(logl, aic, aicc, ak.wts)
        row.names(rt) <- c("GRW", "URW", "Stasis")
        print(rt)
        cat("\n\nParameter estimates: \n")
        print(hats)
    }
    else return(w)
}

