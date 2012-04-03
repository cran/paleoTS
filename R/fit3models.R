fit3models <-
function (y, silent = FALSE, method = c("AD", "Joint"), ...) 
{
    method <- match.arg(method)
    if (method == "AD") {
        m1 <- opt.GRW(y, ...)
        m2 <- opt.URW(y, ...)
        m3 <- opt.Stasis(y, ...)
    }
    else if (method == "Joint") {
        m1 <- opt.joint.GRW(y, ...)
        m2 <- opt.joint.URW(y, ...)
        m3 <- opt.joint.Stasis(y, ...)
    }
    mc <- compareModels(m1, m2, m3, silent = silent)
    invisible(mc)
}
