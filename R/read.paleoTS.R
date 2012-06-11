read.paleoTS <-
function (file = NULL, oldest = "first", reset.time = TRUE, ...) 
{
    if (is.null(file)) {
        ff <- file.choose()
        x <- read.table(ff, ...)
        lab1 <- ff
    }
    else {
        x <- read.table(file = file, ...)
        #lab1 <- paste(getwd(), file)
        lab1<- file
    }
    xr <- as.paleoTS(mm = x[, 2], vv = x[, 3], nn = x[, 1], tt = x[, 
        4], label = lab1, oldest = oldest, reset.time = reset.time)
    return(xr)
}
