## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 5
)

## ------------------------------------------------------------------------
library(paleoTS)
set.seed(10)  # to make example replicatable
x <- sim.GRW(ns = 20, ms = 0.5, vs = 0.1)
plot(x)

## ------------------------------------------------------------------------
print(str(x))

## ------------------------------------------------------------------------
x.sub <- sub.paleoTS(x, ok = 10:15)  # select only populations 10 through 15
plot(x)
plot(x.sub, add = TRUE, col = "red")

## ------------------------------------------------------------------------
 library(mnormt)  # should omit this later
 w.grw <- fitSimple(x, model = "GRW")
 print(w.grw$parameters)  # look at estimated parameters

## ------------------------------------------------------------------------
plot(x, modelFit = w.grw)

## ------------------------------------------------------------------------
 w.urw <- fitSimple(x, model = "URW")
 compareModels(w.grw, w.urw)  # convenient table comparing model support

## ------------------------------------------------------------------------
 w.punc <- fitGpunc(x, ng = 2)  # ng is the number of segments (= number of punctuations + 1)

## ------------------------------------------------------------------------
 plot(x, modelFit = w.punc)

## ------------------------------------------------------------------------
 compareModels(w.grw, w.urw, w.punc)

## ------------------------------------------------------------------------
 fit3models(x)

## ------------------------------------------------------------------------
data(dorsal.spines)
ok1 <- dorsal.spines$nn > 0    # levels without measured fossils
ok2 <- dorsal.spines$tt > 4.4  # levels before the new species invades
ds.sub <- sub.paleoTS(dorsal.spines, ok = ok1 & ok2, reset.time = TRUE)  # subsample 
ds.sub.pool <- pool.var(ds.sub, minN = 5, ret.paleoTS = TRUE)  # replace some pooled variance
w.ou <- fitSimple(ds.sub.pool, pool = FALSE, model = "OU")
plot(ds.sub.pool, modelFit = w.ou)

## ------------------------------------------------------------------------
set.seed(90)
y <- sim.GRW(ns = 40, ms = 0.2, vs = 0.1, vp = 4)  # high vp gives broader error bars
plot(y)

fit3models(y, method = "Joint")  # GRW clearly wins
fit3models(y, method = "AD")     # GRW only barely beats URW


