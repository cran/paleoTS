as.paleoTS <-
function (mm, vv, nn, tt, MM = NULL, genpars = NULL, label = NULL, 
    start.age = NULL, oldest = c("first", "last"), reset.time = TRUE) 
{
    oldest<- match.arg(oldest)
    y <- list(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, genpars = genpars, 
        label = label, start.age = start.age)
    if (oldest == "last") {
        oo <- length(y$mm):1
        y$mm <- y$mm[oo]
        y$vv <- y$vv[oo]
        y$nn <- y$nn[oo]
        y$tt <- y$tt[oo]
    }
    if(y$tt[1] > y$tt[2])	timeDir<- "decreasing"
    else 					timeDir<- "increasing"
    
    if (reset.time) {
        if (y$tt[1] != 0) {
            sa <- y$tt[1]
            if (!is.null(y$start.age) && sa != y$start.age) 
                stop("Age of first sample does not match start.age argument")
            if(timeDir == "decreasing")	y$tt <- sa - y$tt   # decreasing ages (e.g., Ma)  
            else					    y$tt <- y$tt - sa   # increasing ages (elapsed time)
            y$start.age<- sa	
        }
    }
    y$timeDir<- timeDir
    class(y) <- "paleoTS"
    return(y)
}
