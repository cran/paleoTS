`sim.OU` <-
function (ns=20, anc=0, theta=10, alpha=0.3, vs=0.1, vp=1, nn=rep(20, ns), tt=1:ns)
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

