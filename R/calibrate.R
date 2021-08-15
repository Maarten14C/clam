# removed calibrate(), .hpd() and calBP.14C() as these are now provided by the IntCal R package. 
# kept student.t() and .caldist() Aug 2021

# See Christen and Perez 2009, Radiocarbon 51:1047-1059. Instead of assuming the standard Gaussian model (default in clam), a student t distribution can be used with two parameters. Christen and Perez 2009 suggest t.a = 3 and t.b = 4; this can be put as clam( calibt=c(3,4) )
.calibt <- function(t.a, t.b, f.cage, f.error, f.mu, f.sigma) # removed theta as par
  (t.b + ((f.cage-f.mu)^2) / (2*(f.sigma^2 + f.error^2))) ^ (-1*(t.a+0.5))

#(t.b + ((y-x)^2) / (2*(error^2))) ^ (-1*(t.a+0.5))

#' @name student.t 
#' @title Comparison dates calibrated using both the student-t distribution and the the normal distribution.
#' @description Visualise how a date calibrates using the student-t distribution and the the normal distribution.
#' @details Radiocarbon and other dates are usually modelled using the normal distribution (red curve). The student-t approach (grey distribution) however allows for wider tails and thus tends to better accommodate outlying dates. This distribution requires two parameters, called 'a' and 'b'.
#' @param y The reported mean of the date.
#' @param error The reported error of the date.
#' @param t.a Value for the student-t parameter \code{a}.
#' @param t.b Value for the student-t parameter \code{b}.
#' @param postbomb Which postbomb curve to use for negative 14C dates
#' @param cc calibration curve for C14 dates (1, 2 or 3).
#' @param cc1 For northern hemisphere terrestrial C14 dates.
#' @param cc2 For marine C14 dates.
#' @param cc3 For southern hemisphere C14 dates.
#' @param cc4 A custom calibration curve
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are allocated. By default \code{ccdir=""}. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param Cutoff Threshold above which calibrated probabilities are plotted
#' @param times 8 by default.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @author Maarten Blaauw
#' @examples 
#' student.t() 
#' 
#' @export
student.t <- function(y=2450, error=50, t.a=3, t.b=4, cc=1, postbomb=NULL, cc1="IntCal20", cc2="Marine20", cc3="SHCal20", cc4="mixed", ccdir="",Cutoff=1e-5, times=8, rule=1) {
  ccdir <-.validateDirectoryName(ccdir)
  # set the calibration curve
  if(ccdir=="")
    ccdir = paste(system.file("extdata", package="IntCal"), "/", sep="")

  if(cc == 0) {
    x <- seq(y-(times*error), y+(times*error), length=500)
    norm.cal <- dnorm(x, y, error)
    norm.cal <- norm.cal/sum(norm.cal)
    t.cal <- (t.b + ((y-x)^2) / (2*(error^2))) ^ (-1*(t.a+0.5))
    t.cal <- t.cal/sum(t.cal)
    t.cal <- cbind(c(min(x), x, max(x)), c(0, t.cal, 0))
    plot(x, norm.cal, type="l", xlab="cal BP", xlim=range(x)[2:1], ylab="", ylim=c(0, max(t.cal[,2], norm.cal)), col=2, lwd=1.5)
    polygon(t.cal, col=rgb(0,0,0,.25), border=rgb(0,0,0,.5))
    legend("topleft", "Gaussian", text.col=2, bty="n")
    legend("topright", paste("student-t (a=", t.a, ", b=", t.b, ")", sep=""), bty="n", text.col=grey(.4))
  } else {
    if(cc1=="IntCal20") cc1 <- read.table(paste0(ccdir, "3Col_intcal20.14C")) else
      cc1 <- read.csv(paste(ccdir, cc1,  sep=""))[,1:3]
    if(cc2=="Marine20") cc2 <- read.table(paste0(ccdir, "3Col_marine20.14C")) else
      cc2 <- read.csv(paste(ccdir, cc2,  sep=""))[,1:3]
    if(cc3=="SHCal20") cc3 <- read.table(paste0(ccdir, "3Col_shcal20.14C")) else
      cc3 <- read.table(paste(ccdir, cc3,  sep=""))[,1:3]
    if(cc4=="mixed") {
      fl <- paste0(ccdir, "mixed.14C")
      if(file.exists(fl))
        cc4 <- read.table(paste0(ccdir, "mixed.14C")) 
      } else
        cc4 <- read.table(paste0(ccdir, cc4))[,1:3]
    if(cc==1) cc <- cc1 else if(cc==2) cc <- cc2 else if(cc==3) cc <- cc3 else cc <- cc4

    if (y < 0)
      if (length(postbomb) == 0) 
        stop("Warning, negative ages require a postbomb curve. Provide value for postbomb")
    else {
      if(postbomb==1) bomb <- read.table(system.file("extdata","postbomb_NH1.14C", package="IntCal"))[,1:3] else
        if(postbomb==2) bomb <- read.table(system.file("extdata","postbomb_NH2.14C", package="IntCal"))[,1:3] else
          if(postbomb==3) bomb <- read.table(system.file("extdata","postbomb_NH3.14C", package="IntCal"))[,1:3] else
            if(postbomb==4) bomb <- read.table(system.file("extdata","postbomb_SH1-2.14C", package="IntCal"))[,1:3] else
              if(postbomb==5) bomb <- read.table(system.file("extdata","postbomb_SH3.14C", package="IntCal"))[,1:3] else
                stop("Warning, cannot find postbomb curve #", postbomb, " (use values of 1 to 5 only)")
              
      bomb.x <- seq(max(bomb[, 1]), min(bomb[, 1]), length = 500)
      bomb.y <- approx(bomb[, 1], bomb[, 2], bomb.x, rule=rule)$y
      bomb.z <- approx(bomb[, 1], bomb[, 3], bomb.x, rule=rule)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level = 0)
      #if (info$postbomb < 4) #JEV warning
      if (postbomb < 4) 
        cc <- rbind(bomb, cc1, deparse.level = 0)
      else cc <- rbind(bomb, cc3, deparse.level = 0)
    }
    
    norm.cal <- dnorm(cc[, 2], y, sqrt(cc[, 3]^2 + error^2))
    norm.cal <- cbind(cc[, 1], norm.cal/sum(norm.cal))
    acc <- which(norm.cal[, 2] >= Cutoff)
    if(y < 0) acc <- 1:200 # feo pero funciona
    norm.cal <- norm.cal[acc, ]
    
    t.cal <- (t.b + ((y - cc[, 2])^2)/(2 * (cc[, 3]^2 + error^2)))^
      (-1 * (t.a + 0.5))
    t.cal <- cbind(cc[, 1], t.cal/sum(t.cal))
    acc <- which(t.cal[, 2] >= Cutoff)
    if(y < 0) acc <- 1:200 # feo pero funciona
    t.cal <- t.cal[acc, ]
    
    plot(norm.cal, type = "l", xlab = "cal BP", xlim = range(c(t.cal[,1], norm.cal[, 1]))[2:1], ylab = "",
	  ylim = c(0, max(t.cal[,2], norm.cal[, 2])), col = 2, lwd = 1.5)
    polygon(t.cal, col = rgb(0, 0, 0, 0.25), border = rgb(0, 0, 0, 0.5))
    legend("topleft", "Gaussian", text.col = 2, bty = "n")
    legend("topright", paste("student-t (a=", t.a, ", b=", t.b, ")", sep = ""), bty = "n", text.col = grey(0.4))
  }     

}



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


