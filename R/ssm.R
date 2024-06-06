# State-space model functions #


#### model-fitting functions  ####






#' Fit evolutionary models using state-space models (SSM)
#'
#' @param y a \code{paleoTS} object
#' @param z a covariate vector, used only for the covOU models
#' @param gg a grouping vector, used only for the URWshift and covOU_vshift models
#' @param pool if \code{TRUE}, sample variances are substituted with their pooled estimate
#' @param cl optional control list, passed to \code{optim()}
#' @param hess if \code{TRUE}, return standard errors of parameter estimates from the
#' hessian matrix
#'
#' @return a \code{paleoTSfit} object with the model fitting results
#' @export
#' @details These functions use a state space model formulation to compute likelihoods and fit models.
#' @note  For GRW, URW, Stasis, StrictStasis, ACDC and OU models, it will likely be easier to use the convenience function \code{fitSimple} with argument \code{method = "SSM"}.\cr
#' The grouping vector, \code{gg}, is a vector of length equal to the number of samples. It has one element for each sample and takes
#' integer value from 1 to the number of sample groups separated by shiftpoints. See the example below.
#' @seealso \code{\link{fitSimple}}, \code{\link{opt.joint.GRW}}
#' @details
#' Functions to fit the OU covariate tracking models (\code{covOU}, \code{covOU_vshift}) require a covariate argument, \code{z}.
#' At present, only the OU covariate tracking with a shift in the step variance (\code{covOU_vshift}) requires the grouping vector argument (\code{gg}).
#'
#' @examples
#' y <- sim.GRW(ns = 30, vs = 2)
#' w1 <- opt.ssm.URW(y)
#' gg <- rep(1:2, each = 15) # shift occurs immediately after sample 15
#' w2 <- opt.ssm.URWshift(y, gg = gg) # test model in which the step variance shifts
#' compareModels(w1, w2)
opt.ssm.GRW <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates
  p0 <- c(y$mm[1], mle.GRW(y))
  if(p0[3] <= 0) p0[3] <- 1e-5
  names(p0) <- c("anc", "mstep", "vstep")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)

  # do optimization
  fit <- optim(p0, fn = logL.ssm.GRW, method = "L-BFGS-B", control = cl, lower = c(NA, NA, 1e-6), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "GRW", method = "SSM", K = 3,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.GRW")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.GRW(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.URW <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates
  p0 <- c(y$mm[1], mle.URW(y))
  if(p0[2] <= 0) p0[2] <- 1e-5
  names(p0) <- c("anc", "vstep")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)

  # do optimization
  fit <- optim(p0, fn = logL.ssm.URW, method = "L-BFGS-B", control = cl, lower = c(NA, 1e-6), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "URW", method = "SSM", K = 2,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.URW")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.URW(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.Stasis <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates
  p0 <- mle.Stasis(y)
  if(p0[2] <= 0) p0[2] <- 1e-5
  names(p0) <- c("theta", "omega")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)


  # do optimization
  fit <- optim(p0, fn = logL.ssm.Stasis, method = "L-BFGS-B", control = cl, lower = c(NA, 1e-6), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "Stasis", method = "SSM", K = 2,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.Stasis")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.Stasis(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.StrictStasis <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates
  p0 <- mean(y$mm)
  names(p0) <- "theta"
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)


  # do optimization
  fit <- optim(p0, fn = logL.ssm.StrictStasis, method = "L-BFGS-B", control = cl, lower = c(NA), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "StrictStasis", method = "SSM", K = 1,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.StrictStasis")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.StrictStasis(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.OU <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates; start with alpha giving half-life equal to half the sequence
  p0 <- c(y$mm[1], mle.URW(y), mean(y$mm), log(2)/(0.5*max(y$tt)) )
  if(p0[2] <= 0) p0[2] <- 1e-3
  names(p0) <- c("anc", "vstep", "theta", "alpha")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)


  # do optimization
  fit <- optim(p0, fn = logL.ssm.OU, method = "L-BFGS-B", control = cl, lower = c(NA, 1e-6, NA, 1e-6), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "OU", method = "SSM", K = 4,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.OU")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.OU(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.ACDC <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){
  ## check if pooled, make start at tt=0
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
  if(y$tt[1] != 0)	stop("Initial time should be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")

  # get initial estimates; start at URW [r = 0]
  p0 <- c(y$mm[1], mle.URW(y), 0)
  if(p0[2] <= 0) p0[2] <- 1e-3
  names(p0) <- c("anc", "vstep0", "r")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)


  # do optimization
  fit <- optim(p0, fn = logL.ssm.ACDC, method = "L-BFGS-B", control = cl, lower = c(NA, 1e-6, NA), hessian = hess, y = y)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "ACDC", method = "SSM", K = 3,
                     n = length(y$mm), se = se, convergence = fit$convergence, logLFunction = "logL.ssm.ACDC")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.ACDC(p = w$parameters, y, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.covOU <- function(y, z, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){

  # handle covariate: optimum position is theta = b0 + b1*z
  ns <- length(y$mm)
  if (length(z) != ns)
    stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n")
  zc <- z - mean(z)
  cat("Note: mean centering the covariate (z) to facilitate model-fitting.\n")

  # get initial estimates; start with alpha giving half-life equal 5 avg steps
  wr <- stats::lm(y$mm ~ zc)
  #p0 <- c(y$mm[1], mle.URW(y), log(2)/(5*mean(diff(y$tt))), coef(wr)[1], coef(wr)[2])
  p0 <- c(y$mm[1], mle.URW(y), log(2)/(5*mean(diff(y$tt))), mean(y$mm), 0)
  if(p0[3] <= 1e-5) p0[3] <- 1e-5 # fine for dt = 1, but probably want a term for avg dt
  if(p0[2] <= 1e-5) p0[2] <- 1e-5
  names(p0) <- c("anc", "vstep", "alpha", "b0", "b1")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0, z)

  # do optimization; if fails, parscale is probably the problem
  fit <- optim(p0, fn = logL.ssm.covOU, method = "L-BFGS-B", lower = c(NA, 1e-6, 1e-6, NA, NA),
               control=cl, hessian = hess, y = y, z = zc)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL

  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "covOU",
                     method = "SSM", K = length(fit$par), n = length(y$mm), se = se,
                     convergence = fit$convergence, logLFunction = "logL.ssm.covOU")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.covOU(p = w$parameters, y, zc, logL.only = FALSE)

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.URWshift <- function(y, gg, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){

  nseg <- max(gg) # number of segments
  if(nseg < 2) stop("This model requires at least 2 segments.")

  # get initial estimates
  vs0 <- mle.URW(y)
  if(vs0 <= 0) vs0 <- 1e-3
  p0 <- c(y$mm[1], rep(vs0, nseg))
  names(p0) <- c("anc", paste0("vstep", 1:nseg))
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)

  # do optimization
  ll <- c(NA, rep(1e-6, nseg))
  fit <- optim(p0, fn = logL.ssm.URWshift, method = "L-BFGS-B", lower = ll, control = cl,
               hessian = hess, y = y, gg = gg)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  mn <- paste0("URW-shift", nseg-1)
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = mn,
                     method = "SSM", K = length(p0), n = length(y$mm), se = se,
                     convergence = fit$convergence, logLFunction = "logL.ssm.URWshift")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.URWshift(p = w$parameters, y, gg, logL.only = FALSE)
  w$gg <- gg
  #w$p0 <- p0

  return(w)
}

#' @rdname opt.ssm.GRW
#' @export
opt.ssm.covOU_vshift <- function(y, z, gg, pool = TRUE, cl = list(fnscale = -1), hess = FALSE){

  # handle covariate: version where theta = b0 + b1*z
  # also possible: d theta = b1*dz
  ns <- length(y$mm)
  if (length(z) != ns)
    stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n")
  zc <- z - mean(z)
  cat("Note: mean centering the covariate (z) to facilitate model-fitting.\n")


  # get initial estimates; start with alpha giving half-life equal to half the sequence
  nseg <- max(gg) # number of segments
  if(nseg < 2) stop("This model requires at least 2 segments.")

  # initially, take all initial vs as equal =
  wr <- stats::lm(y$mm ~ zc)
  vs0 <- mle.URW(y)
  if(vs0 <= 0) vs0 <- 1e-3
  p0 <- c(y$mm[1], rep(vs0, nseg), log(2)/(5*mean(diff(y$tt))), stats::coef(wr)[1], stats::coef(wr)[2])
  names(p0) <- c("anc", paste0("vstep", 1:nseg), "alpha", "b0", "b1")
  if(is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0, z)

  # do optimization
  fit <- optim(p0, fn = logL.ssm.covOU_vshift, method = "L-BFGS-B",
               lower = c(NA, rep(1e-6, nseg), 1e-6, NA, NA),
               control=cl, hessian = hess, y = y, z = zc, gg = gg)

  # gather information to return
  if (hess){
    se <- sqrt(diag(-1 * solve(fit$hessian)))
  } else se <- NULL
  w <- as.paleoTSfit(logL = fit$value, parameters = fit$par, modelName = "covOU_vshift",
                     method = "SSM", K = length(fit$par), n = length(y$mm), se = se,
                     convergence = fit$convergence, logLFunction = "logL.ssm.covOU_vshift")

  # add Kalman filter, computed at the MLE, to the paleoTSfit object returned
  w$kf <- logL.ssm.covOU_vshift(p = w$parameters, y, zc, gg, logL.only = FALSE)
  w$convergence <- fit$convergence
  #w$counts <- fit$counts
  #w$p0 <- p0
  w$gg <- gg

  return(w)
}



#### unexported logL functions ####

# URW log-likelihood function via Kalman filter, not exported
logL.ssm.URW <- function(p, y, logL.only = TRUE){
  anc <- p[1]
  vs <- p[2]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(1, dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array(vs*dt, dim = c(1,1,n))

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = NULL, Gam = NULL, Qtv, Rtv, input = NULL)
  if(logL.only) return(kf$like) else return(kf)
}

# GRW log-likelihood function via Kalman filter, not exported
logL.ssm.GRW <- function(p, y, logL.only = TRUE){
  anc <- p[1]
  ms <- p[2]
  vs <- p[3]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(1, dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array(vs*dt, dim = c(1,1,n))
  input <- ms*dt

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = 1, Gam = NULL, Qtv, Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}

logL.ssm.Stasis <- function(p, y, logL.only = TRUE){
  theta <- p[1]
  omega <- p[2]
  n <- length(y$mm) # number of pop's in time-series
  #dt <- c(0, diff(y$tt))

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(0, dim = c(1,1,n))
  mu0 <- y$mm[1]
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array(omega, dim = c(1,1,n))
  input <- rep(1, n)

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = theta, Gam = NULL, Qtv, Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}

logL.ssm.StrictStasis <- function(p, y, logL.only = TRUE){
  theta <- p[1]
  n <- length(y$mm) # number of pop's in time-series

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(0, dim = c(1,1,n))
  mu0 <- y$mm[1]
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array(0, dim = c(1,1,n))
  input <- rep(1, n)

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = theta, Gam = NULL, Qtv, Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}


logL.ssm.OU <- function(p, y, logL.only = TRUE){
  anc <- p[1]
  theta <- p[3]  # changed order to match Joint function
  vs <- p[2]
  alpha <- p[4]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(exp(-alpha*dt), dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array( (vs/(2*alpha)) * (1 - exp(-2*alpha*dt)), dim = c(1,1,n))
  input <- theta * (1 - exp(-alpha*dt))

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = 1, Gam = NULL, Qtv, Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}

logL.ssm.ACDC <- function(p, y, logL.only = TRUE){
  anc <- p[1]
  vs <- p[2]
  r <- p[3]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))
  si <- y$tt - y$tt[1] # elapsed time, in case tt[1] != 0

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(1, dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))

  #state variance, need to handle when r = 0 to avoid NaN
  if(r == 0) Qtv <- array(vs*dt, dim = c(1,1,n)) else{
    V <- vs*(exp(r*si) - 1)/r
    Qtv <- array(c(0, diff(V)), dim = c(1,1,n))
  }

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = NULL, Gam = NULL, Qtv, Rtv, input = NULL)
  if(logL.only) return(kf$like) else return(kf)
}


# log-likelihood function for OU-cov-tracking
logL.ssm.covOU <- function(p, y, z, logL.only = TRUE){
  anc <- p[1]
  vs <- p[2]
  alpha <- p[3]
  b0 <- p[4]
  b1 <- p[5]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(exp(-alpha*dt), dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array( (vs/(2*alpha)) * (1 - exp(-2*alpha*dt)), dim = c(1,1,n))
  input <- (b0 + b1*z) * (1 - exp(-alpha*dt)) # theta replaced by f(z)

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}

# log-likelihood function for URW with shift in vstep
logL.ssm.URWshift <- function(p, y, gg, logL.only = TRUE){
  nseg <- max(gg)
  anc <- p[1]
  vsv <- p[2:(nseg+1)]
  vss <- vsv[gg]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))
  if(length(vss)==2) print(vss)

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(1, dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array(vss*dt, dim = c(1,1,n))

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = NULL, Gam = NULL, Qtv, Rtv, input = NULL)
  if(logL.only) return(kf$like) else return(kf)
}

# log-likelihood for covOU with shift in vstep
logL.ssm.covOU_vshift <- function(p, y, z, gg, logL.only = TRUE){
  nseg <- max(gg)
  anc <- p[1]
  mxv <- nseg + 1
  vsv <- p[2:mxv]
  alpha <- p[mxv + 1]
  b0 <- p[mxv + 2]
  b1 <- p[mxv + 3]
  n <- length(y$mm) # number of pop's in time-series
  dt <- c(0, diff(y$tt))
  vss <- vsv[gg] # vector of vsteps for each of n samples

  # set-up Kfilter elements
  Atv <- array(1, dim = c(1,1,n))
  Phitv <- array(exp(-alpha*dt), dim = c(1,1,n))
  mu0 <- anc
  Sigma0 <- y$vv[1]/y$nn[1]
  Rtv <- array(y$vv/y$nn, dim = c(1,1,n))
  Qtv <- array( (vss/(2*alpha)) * (1 - exp(-2*alpha*dt)), dim = c(1,1,n))
  input <- (b0 + b1*z) * (1 - exp(-alpha*dt)) # theta replaced by f(z)

  kf <- Kfiltertv(n, y$mm, Atv, mu0, Sigma0, Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input)
  if(logL.only) return(kf$like) else return(kf)
}






#' Time-varying Kalman filter calculations
#'
#' @param num the number of samples in the time-series
#' @param y values of the time-series
#' @param Atv q x p x n observation array
#' @param mu0 p x 1 vector setting the mean of the system at time zero
#' @param Sigma0 p x p variance matrix of the system at time zero
#' @param Phitv p x p x n array reflecting autoregression of the state variables
#' @param Ups p x r matrix with the coefficients/parameters relating the inputs to the system equation
#' @param Gam q x r matrix with the coefficients/parameters relating the inputs to the observation equation
#' @param Qtv p x p x n array of system stochastic  variance; user needs to ensure positive definite
#' @param Rtv q x q x n array observation stochastic variance; user needs to ensure positive definite
#' @param input n x r array of the exogenous variables/covariates
#'
#' @return A list of the following elements:
#' \itemize{
#' \item \code{xp}   one-step-ahead prediction of the state
#' \item \code{Pp}   mean square prediction error
#' \item \code{xf}   filter value of the state
#' \item \code{Pf}   mean square filter error
#' \item\code{like}   log-likelihood
#' \item\code{innov}   innovation series
#' \item\code{sig}   innovation covariances
#' \item\code{Kn}   last value of the gain, needed for smoothing
#' }
#' @details
#' For the dimensions of the argument arrays, \code{n} is the length of the
#' time-series, \code{q} is the dimension of the observation variable(s),
#' \code{p} is the dimension of the state variable(s), and \code{r} isthe
#' dimension of the input variable(s). \cr\cr This function is based on the
#' \code{Kfilter} function of the \emph{astsa} package, modified modified to
#' allow for time-varying terms for the Kalman filter.  This modification
#' facilitates fitting a broader array of models and handling non-uniform
#' temporal spacing of samples.  See the documentation for that function, and
#' the reference below for additional information.
#'
#' @note This function is used in the internal SSM log-likelihood functions for the models.  The user will not need to use this
#' they create their own model-fitting functions.
#' @author John Fricks (jfricks@asu.edu)
#' @references Shumway, R. H., and D. S. Stoffer. 2017. Time Series Analysis and its Applications (4th Ed.) Springer International.
#'
#' @export
#' @examples
#' y <- sim.GRW(ms = 0, vs = 1, vp = 0)
#' n <- length(y)
#' kf <- Kfiltertv(n ,y = y$mm, Atv = array(1, dim = c(1,1,n)), mu0 = y$mm[1],
#'                 Sigma0 = y$vv[1]/y$nn[1], Phitv = array(1, dim = c(1,1,n)),
#'                 Ups = NULL, Gam = NULL, Qtv = array(1, dim = c(1,1,n)),
#'                 Rtv = array(0, dim = c(1,1,n)), input = NULL)
#'
Kfiltertv <-function(num,y,Atv,mu0,Sigma0,Phitv,Ups,Gam,Qtv,Rtv,input){
  if (is.null(input)){  # no input
    Q = as.matrix(Qtv[,,1])
    R=as.matrix(Rtv[,,1])
    Phi    = as.matrix(Phitv[,,1])
    pdim   = nrow(Phi)
    mu0    = matrix(mu0, nrow=pdim, ncol=1)
    Sigma0 = matrix(Sigma0, nrow=pdim, ncol=pdim)
    y      = as.matrix(y)
    qdim   = ncol(y)
    xp    = array(NA, dim=c(pdim,1,num))         # xp = x_t^{t-1}
    Pp    = array(NA, dim=c(pdim,pdim,num))      # Pp = P_t^{t-1}
    xf    = array(NA, dim=c(pdim,1,num))         # xf = x_t^t
    Pf    = array(NA, dim=c(pdim,pdim,num))      # Pf = x_t^t
    innov = array(NA, dim=c(qdim,1,num))         # innovations
    sig   = array(NA, dim=c(qdim,qdim,num))      # innov var-cov matrix
    # initialize
    xp[,,1]  = Phi%*%mu0
    Pp[,,1]  = Phi%*%Sigma0%*%t(Phi)+Q
    B       = matrix(Atv[,,1], nrow=qdim, ncol=pdim)
    sigtemp  = B%*%Pp[,,1]%*%t(B) + R
    sig[,,1] = (t(sigtemp)+sigtemp)/2    # innov var - make sure symmetric
    siginv   = solve(sig[,,1])
    K        = Pp[,,1]%*%t(B)%*%siginv
    innov[,,1] = y[1,] - B%*%xp[,,1]
    xf[,,1]  = xp[,,1] + K%*%innov[,,1]
    Pf[,,1]  = Pp[,,1] - K%*%B%*%Pp[,,1]
    sigmat   = matrix(sig[,,1], nrow=qdim, ncol=qdim)
    like     = log(det(sigmat)) + t(innov[,,1])%*%siginv%*%innov[,,1]   # -log(likelihood)

    # start filter iterations

    for (i in 2:num){
      if (num < 2) break
      Q = as.matrix(Qtv[,,i])
      R=as.matrix(Rtv[,,i])
      Phi    = as.matrix(Phitv[,,i])
      xp[,,i]   = Phi%*%xf[,,i-1]
      Pp[,,i]   = Phi%*%Pf[,,i-1]%*%t(Phi)+Q
      B        = matrix(Atv[,,i], nrow=qdim, ncol=pdim)
      sigma     = B%*%Pp[,,i]%*%t(B) + R
      sig[,,i]  = (t(sigma)+sigma)/2     # make sure sig is symmetric
      siginv    = solve(sig[,,i])          # now siginv is sig[[i]]^{-1}
      K         = Pp[,,i]%*%t(B)%*%siginv
      innov[,,i] = y[i,] - B%*%xp[,,i]
      xf[,,i]   = xp[,,i] + K%*%innov[,,i]
      Pf[,,i]   = Pp[,,i] - K%*%B%*%Pp[,,i]
      sigmat    = matrix(sig[,,i], nrow=qdim, ncol=qdim)
      like      = like + log(det(sigmat)) + t(innov[,,i])%*%siginv%*%innov[,,i]
    }  # end no input
    #
  } else {  # start with input
    #
    Q = as.matrix(Qtv[,,1])
    R=as.matrix(Rtv[,,1])
    Phi    = as.matrix(Phitv[,,1])
    pdim   = nrow(Phi)
    mu0    = matrix(mu0, nrow=pdim, ncol=1)
    Sigma0 = matrix(Sigma0, nrow=pdim, ncol=pdim)
    y      = as.matrix(y)
    qdim   = ncol(y)
    rdim   = ncol(as.matrix(input))
    input  = matrix(input, nrow=num, ncol=rdim)
    if (is.null(Ups)) Ups = matrix(0, nrow=pdim, ncol=rdim)
    if (is.null(Gam)) Gam = matrix(0, nrow=qdim, ncol=rdim)
    Ups = as.matrix(Ups)
    Gam = as.matrix(Gam)
    xp    = array(NA, dim=c(pdim,1,num))         # xp = x_t^{t-1}
    Pp    = array(NA, dim=c(pdim,pdim,num))      # Pp = P_t^{t-1}
    xf    = array(NA, dim=c(pdim,1,num))         # xf = x_t^t
    Pf    = array(NA, dim=c(pdim,pdim,num))      # Pf = x_t^t
    innov = array(NA, dim=c(qdim,1,num))         # innovations
    sig   = array(NA, dim=c(qdim,qdim,num))      # innov var-cov matrix
    # initialize
    xp[,,1]  = Phi%*%mu0 + Ups%*%input[1,]
    Pp[,,1]  = Phi%*%Sigma0%*%t(Phi)+Q
    B       = matrix(Atv[,,1], nrow=qdim, ncol=pdim)
    sigtemp  = B%*%Pp[,,1]%*%t(B) + R
    sig[,,1] = (t(sigtemp)+sigtemp)/2    # innov var - make sure symmetric
    siginv   = solve(sig[,,1])
    K        = Pp[,,1]%*%t(B)%*%siginv
    innov[,,1] = y[1,] - B%*%xp[,,1] - Gam%*%input[1,]
    xf[,,1]  = xp[,,1] + K%*%innov[,,1]
    Pf[,,1]  = Pp[,,1] - K%*%B%*%Pp[,,1]
    sigmat   = matrix(sig[,,1], nrow=qdim, ncol=qdim)
    like     = log(det(sigmat)) + t(innov[,,1])%*%siginv%*%innov[,,1]   # -log(likelihood)


    # start filter iterations

    for (i in 2:num){
      if (num < 2) break
      Q = as.matrix(Qtv[,,i])
      R=as.matrix(Rtv[,,i])
      Phi    = as.matrix(Phitv[,,i])
      xp[,,i]   = Phi%*%xf[,,i-1] + Ups%*%input[i,]
      Pp[,,i]   = Phi%*%Pf[,,i-1]%*%t(Phi)+Q
      B        = matrix(Atv[,,i], nrow=qdim, ncol=pdim)
      sigma     = B%*%Pp[,,i]%*%t(B) + R
      sig[,,i]  = (t(sigma)+sigma)/2     # make sure sig is symmetric
      siginv    = solve(sig[,,i])          # now siginv is sig[[i]]^{-1}
      K         = Pp[,,i]%*%t(B)%*%siginv
      innov[,,i] = y[i,] - B%*%xp[,,i] - Gam%*%input[i,]
      xf[,,i]   = xp[,,i] + K%*%innov[,,i]
      Pf[,,i]   = Pp[,,i] - K%*%B%*%Pp[,,i]
      sigmat    = matrix(sig[,,i], nrow=qdim, ncol=qdim)
      like      = like + log(det(sigmat)) + t(innov[,,i])%*%siginv%*%innov[,,i]
    }
  } # end with input
  like = -0.5*like # note: changed sign here to return logL, not -logL, to keep consistent with other paleoTS functions
  list(xp=xp,Pp=Pp,xf=xf,Pf=Pf,like=like,innov=innov,sig=sig,Kn=K)
}



#' Compute and (optionally) plot residuals from SSM model fit
#'
#' @param y a \code{paleoTS} object
#' @param w a \code{paleoTSfit} object, where the fitting functions uses SSM
#' @param show.plot logical, if \code{TRUE}, several diagnostic plots are shown
#' @param resid.type residual type, either "standardized" or "unstandardized"
#' @details
#' It is recommended that \code{resid.type} be set to the default, "standardized", which will scale residuals by their expected standard deviation
#'
#' @return a vector of residuals, returned invisibly
#' @export
#'
#' @examples
#' y <- sim.GRW(ns = 50, ms = 0.2)
#' w <- fitSimple(y, model = "URW", method = "SSM")  # wrong model
#' checkSSMresiduals(y, w, show.plot = TRUE) # positive residuals show model mis-fit
checkSSMresiduals <- function(y, w, show.plot = TRUE, resid.type = c("standardized", "unstandardized")){
  resid.type <- match.arg(resid.type)
  pred <- w$kf$xp
  rr <- y$mm - pred[1,1,]

  if(resid.type == "standardized") rr <- rr / (sqrt(w$kf$Pp[1,1,]))

  # optionally, show diagnostic plots
  if(show.plot){
    graphics::par(ask = TRUE)
    yl <- paste(resid.type, "residuals", sep = " ")
    plot(pred, rr, ylab = yl, xlab = "predicted values")
    mtext("Residuals vs Predicted")
    abline(h=0, lty=3)
    ll <- stats::lowess(pred, rr)
    lines(ll, col = "red")

    stats::qqnorm(rr)
    mtext(yl)
    stats::qqline(rr)

    # now plot residuals in time-order
    if(y$timeDir == "decreasing") real_tt <- y$start.age - y$tt else real_tt <- y$tt

    plot(real_tt, rr, xlab = "time", ylab = "residuals", xlim = rev(range(real_tt)))
    mtext("Residuals over Time")
    abline(h = 0, lty = 3)
    ll <- stats::lowess(real_tt, rr)
    lines(ll, col = "red")

    graphics::par(ask = FALSE)
  }

  invisible(rr)
}

