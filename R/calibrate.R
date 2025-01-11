# functions such as calibrate(), .hpd() and calBP.14C() are now provided by the rice R package.
# kept calib.t() and .caldist() Aug 2021
# removed calib.t() Jan 2025, since it is in `rice` already


# See Christen and Perez 2009, Radiocarbon 51:1047-1059. Instead of assuming the standard Gaussian model (default in clam), a student t distribution can be used with two parameters. Christen and Perez 2009 suggest t.a = 3 and t.b = 4; this can be put as clam( calibt=c(3,4) )
.calibt <- function(t.a, t.b, f.cage, f.error, f.mu, f.sigma) # removed theta as par
  (t.b + ((f.cage-f.mu)^2) / (2*(f.sigma^2 + f.error^2))) ^ (-1*(t.a+0.5))


# find the calibrated distributions of 14C dates
.caldist <- function(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD, normalise=FALSE, rule=rule) {
    if(f.cage > 1) {
      if(f.cage > 1) yrsteps <- min(yrsteps, .1)
        pb <- theta[which(f.mu > 1)]
      if(length(pb)==0)
        stop("help, something exploded with a postbomb date")
      x <- approx(theta, f.mu, seq(min(pb), max(pb), by=yrsteps), rule=2) # Nov 2020 added rule=1
      xsd <- approx(theta, f.sigma, x$x, rule=2)$y # Nov 2020 added rule=1
      theta <- c(x$x, theta[which(f.mu <= 0)])
      f.mu <- c(x$y, f.mu[which(f.mu <= 0)])
      f.sigma <- c(xsd, f.sigma[which(f.mu <= 0)])
      threshold <- 0
    }

    # calibrate; find how far f.cage (measurement) is from f.mu (calibration curve)
    if(length(calibt) < 2)
      cal <- cbind(theta, dnorm(f.mu, f.cage, sqrt(f.error^2+f.sigma^2))) else
        cal <- cbind(theta, .calibt(calibt[1], calibt[2], f.cage, f.error, f.mu, f.sigma))

    # interpolate and normalise calibrated distribution to 1
    cal <- cal[min(which(cal[,2] > 0)):max(which(cal[,2] > 0)),] # remove unnecessary data
    cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=yrsteps), rule=rule)
    cal <- cbind(cal$x, cal$y/sum(cal$y))
    if(BCAD && (0 %in% cal[,1]))
      cal <- cal[-which(cal[,1]==0),] # 0 BC/AD does not exist
    # only report those normalised calibrated probabilities beyond a threshold
    cal[cal[,2] > threshold,]
  }


