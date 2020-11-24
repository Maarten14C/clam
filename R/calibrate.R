
#' @name calibrate 
#' @title Calibrate individual 14C dates.
#' @description Calibrate individual 14C dates, plot them and report calibrated ranges.
#' @details 
#' Type \code{calibrate()} to see how a date of 2450 +- 50 14C BP gets calibrated (the calibration curve happens to show
#' a plateau around this 14C age). To calibrate a different date, provide its reported mean and error (1 
#' standard deviation error as reported by the radiocarbon laboratory) as follows: \code{calibrate(mean, error)},
#' e.g., for a date of 130 +- 20 14C BP, type calibrate\code{(cage=130, error=20)} or, shorter, \code{calibrate(130,20)}. 
#' As this date will fall partly beyond the younger extreme of the calibration curve, a warning will be given
#' (similar warnings will be given for too old dates).
#'   
#' In case the date has a reservoir effect or age offset, e.g. of 100 14C years, provide this as follows: 
#' \code{calibrate(130, 20, reservoir=100)}. If you want to include an uncertainty for this offset, provide this as follows,
#' e.g., for an uncertainty of 50yr, \code{calibrate(130,20,reservoir=c(100, 50))}. 
#' The uncertainty for the age offset will then be added to the error (by taking the square root of the sum 
#' of the squared error and the squared offset uncertainty). If the carbon of your sample has mixed marine/terrestrial sources,
#' instead apply the marine offset using mix.calibrationcurves mix.calibrationcurves, and calibrate the date using that custom-built curve.
#'
#' If you prefer to work with, e.g., 68 \% as opposed to the default 95 \% confidence intervals, 
#' type: \code{calibrate(130, 20, prob=0.68)} or \code{calibrate(130, 20,, 0.68)} (the commas between the brackets indicate the position of the option;
#' the standard deviation is the fourth option of the \code{calibrate} function). Clam calculates the calibrated distribution 
#' for every single calendar year (\code{yrsteps=1}) within a wide range of the 14C date (default but adaptable \code{times=5}
#' standard deviations or 99.999999 \% of its probability distribution). This range can also be adapted by 
#' changing the option expand (default \code{expand=0.1}). Probabilities below a threshold (default \code{threshold=1e-6}) will be neglected.
#' 
#' By default the northern hemisphere terrestrial calibration curve is used (\code{cc=1, cc1="3Col_intcal20.14C"}). 
#' To use alternative curves, use \code{cc=2} (\code{cc2="3Col_marine20.14C"}), \code{cc=3} (\code{cc3="3Col_shcal20.14C"}), 
#' \code{cc=4} (\code{cc4="mixed.14C"}), or change the file names of \code{cc1, cc2, cc3 or cc4}.
#'  
#' Clam works in cal BP (calendar years before AD 1950) by default, but can work with cal BC/AD through the option \code{BCAD=TRUE}. 
#' 
#' By default the Gaussian distribution is used to calibrate dates. For use of the student-t distribution instead, 
#' provide two sensible values, e.g., \code{calibt=c(3,4)}.
#'
#' Calibrated distributions are usually reduced to their 68\% or 95\% calibrated ranges, taking into account the asymmetric 
#' and multi-peaked shape of these distributions. In clam, this is done by calculating the highest posterior density (hpd) ranges: 
#' \itemize{
#' \item i) the probability distribution (see above) is normalised to 100\%
#' \item ii) the calendar years are ranked according to their probabilities
#' \item iii) those calendar ages with a cumulative sum at or above the desired probability threshold (default 95\%) are retained, and 
#' \item iv) the extremes and probabilities of any sub-ranges within these calendar ages are reported. 
#' }
#' Calibrated ranges at 68\% will obviously result in narrower confidence intervals, and a perceived higher precision, than 95\% ranges. However, given the often
#' asymmetric and multi-modal nature of calibrated distributions, the probability that the 'true' calendar date 
#' lies outside the 1 standard deviation hpd ranges is considerable (c. 32\%). Therefore the use of 95\% calibrated ranges is preferable, 
#' and default in clam. The hpd ranges are calculated at yearly resolution by default (\code{hpdsteps=1}).
#'
#' Negative radiocarbon ages are calibrated with postbomb curves, but the user needs to tell clam which curve to use. 
#' For example, to use the first of the three northern hemisphere curves, provide the option \code{postbomb=1}, 
#' while for southern hemisphere samples, use \code{postbomb=4} or \code{postbomb=5}. Default curves can be changed; 
#' currently they are \code{pb1="postbomb_NH1.14C"}, \code{pb2="postbomb_NH2.14C"}, \code{pb3="postbomb_NH3.14C"}, 
#' \code{pb4="postbomb_SH1-2.14C"} and \code{pb5="postbomb_SH3.14C"}; see \url{http://calib.org/CALIBomb/}. 
#' If no \code{postbomb} option is provided 
#' for negative radiocarbon ages, clam will report an error and refuse to calibrate the date. Given the sub-year resolution of postbomb-curves, 
#' hpd ranges are calculated at high resolution by default (\code{pbsteps=0.01}). Choose alternative values with care as 
#' they may cause unexpected results. 
#' 
#' Generally the calculations are removed from memory after calibration; 
#' if you want to have them stored (say for subsequent manipulations), provide the option \code{storedat=TRUE}.
#'
#' A graph of the calibration is produced by default (\code{graph=TRUE}), and it can be adapted in several ways.
#' The limits of the horizontal (calendar scale) and vertical (14C scale) axes are calculated automatically 
#' but can be changed by providing alternative values for the options \code{yrmin, yrmax, minC14} and \code{maxC14}, respectively.
#' The titles of both axis can be changed by providing alternative titles to \code{xlab} and/or \code{ylab}, and 
#' also the top title can be adapted using title. The heights of the distributions of the 14C and calibrated 
#' ages can be set to alternative values using \code{calheight} (default \code{0.3} which plots the distribution up to 30\% of the height of the entire graph).
#' Parameters for white space around the 
#' graph can be changed (default \code{mar=c(3.5, 2, 2, 1}) for spacing below, to the left, above and to the right respectively), 
#' as can the spacing for the axis labels (\code{mgp=c(2,1,0)}). By default, the axes are connected at the lower left, \code{bty="l"}.
#' Check the R documentation of \code{par()} for more options.
#'   
#' The colours of the 14C date, the calibration curve, the entire distributions, as well as of the highest posterior density (\code{hpd}) 
#' ranges, can be changed by providing an alternative colour in \code{date.col}, \code{cc.col}, \code{dist.col}, and/or \code{sd.col}, respectively.
#' The default colours are transparent grey for the dates probability distributions (\code{dist.col=rgb(0,0,0, 0.3)} and \code{sd.col=rgb(0,0,0, 0.5)};
#' change the last value of rgb for different greyscale values), red for the uncalibrated mean and error bars (\code{date.col="red"}), 
#' and transparent green for the calibration curve (\code{cc.col=rgb(0, 0.5, 0, 0.7)}). R's rgb() function expects values between \code{0} and \code{1}
#' for red, green and blue, respectively, followed by a value for the semi-transparency (also between 0 and 1). Some graphic devices 
#' such as postscript are unable to use transparency; in that case provide different colours or leave the fourth value empty.
#' @param cage Mean of the uncalibrated C-14 age.
#' @param error Error of the uncalibrated C-14 age.
#' @param reservoir Reservoir age, or reservoir age and age offset.
#' @param prob Probability confidence intervals (between 0 and 1).
#' @param cc Calibration curve for C-14 dates (1, 2, 3, or 4).
#' @param cc1 For northern hemisphere terrestrial C-14 dates.
#' @param cc2 For marine C-14 dates.
#' @param cc3 For southern hemisphere C-14 dates.
#' @param cc4 For mixed marine/terrestrial C-14 dates.
#' @param ccdir Directory where the calibration curves for C-14 dates \code{cc} are located. By default \code{ccdir=""}. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="/Curves"} to choose sub-folder \code{/Curves}.
#' @param postbomb Calibration curve for postbomb dates.
#' @param pb1 For Northern hemisphere region 1 postbomb C-14 dates.
#' @param pb2 For Northern hemisphere region 2 postbomb C-14 dates.
#' @param pb3 For Northern hemisphere region 3 postbomb C-14 dates.
#' @param pb4 For Southern hemisphere regions 1-2 postbomb C-14 dates.
#' @param pb5 For Southern hemisphere region 3 postbomb C-14 dates.
#' @param yrsteps Temporal resolution at which C-14 ages are calibrated (in calendar years).
#' @param pbsteps Temporal resolution at which postbomb C-14 ages are calibrated (in calendar years).
#' @param hpdsteps Temporal resolution at which highest posterior density ranges are calibrated (in calendar years).
#' @param calibt Calibration based on the student-t distribution. By default, the Gaussian distribution is used (\code{calibt=FALSE}). To use the student-t distribution, provide two parameters such as \code{calibt=c(3,4)}.
#' @param yrmin Minimum of calendar axis (default calculated automatically).
#' @param yrmax Maximum of calendar axis (default calculated automatically).
#' @param minC14 Minimum age of the C-14 age axis (default calculated automatically).
#' @param maxC14 Maximum of the C-14 age axis (default calculated automatically).
#' @param times Half-range of calibration curve used to calibrate dates (multiplication factor for the date's errors).
#' @param calheight Maximum height of the C14 and calibrated distributions (as proportion of the invisible secondary axes).
#' @param expand By which ratio should the calendar axis be expanded to fit the calibrated distribution.
#' @param threshold Below which value should probabilities be excluded from calculations.
#' @param graph Plot a graph of the calibrated date. If set to FALSE, only the hpd ranges will be given.
#' @param storedat Store the dates within the R session after a clam run.
#' @param xlab Label of the horizontal axis. Defaults to the calendar scale, but alternative names can be provided.
#' @param ylab Label of the vertical axis. Defaults to the 14C scale, but alternative names can be provided.
#' @param BCAD Use BC/AD or cal BP scale (default cal BP).
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param bty Draw a box around the graph ("n" for none, and "l", "7", "c", "u", "]" or "o" for correspondingly shaped boxes).
#' @param xaxs Whether or not to extend the limits of the horizontal axis. Defaults to \code{xaxs="i"} which does not extend the limits.
#' @param yaxs Whether or not to extend the limits of the vertical axis. Defaults to \code{yaxs="i"} which does not extend the limits.
#' @param title	Title of the graph. Defaults to the values of the uncalibrated date.
#' @param date.col Colour of the "dot-bar" plot of the C14 date. Defaults to \code{date.col="red"}.
#' @param cc.col Colour of the calibration curve. Defaults to semi-transparent dark green; \code{cc.col=rgb(0,.5,0,0.7)}.
#' @param dist.col Colour of the calibrated distribution.
#' @param sd.col Colour of calibrated range.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @author Maarten Blaauw
#' @return A graph of the raw and calibrated C-14 date, the calibrated ranges and, invisibly, the calibrated ranges and probabilities.
#' @examples 
#' calibrate()
#' calibrate(130, 20)
#' cal <- calibrate(2550, 20, reservoir=100)
#' cal; plot(cal$calib)
#' calibrate(130, 20, prob=0.68)
#' calibrate(cage=130, error=20)
#' calibrate(4450, 40, reservoir=c(100, 50))
#' @export
calibrate <- function(cage=2450, error=50, reservoir=0, prob=0.95, cc=1, cc1="3Col_intcal20.14C", cc2="3Col_marine20.14C", cc3="3Col_shcal20.14C", cc4="mixed.14C", ccdir="", postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C", pb5="postbomb_SH3.14C", yrsteps=1, pbsteps=0.01, hpdsteps=1, calibt=FALSE, yrmin=NULL, yrmax=NULL, minC14=NULL, maxC14=NULL, times=5, calheight=0.3, expand=0.1, threshold=1e-6, storedat=FALSE, graph=TRUE, xlab=NULL, ylab=NULL, BCAD=FALSE, mar=c(3.5,3,2,1), mgp=c(1.7,.8,0), bty="l", xaxs="i", yaxs="i", title=NULL, date.col="red", cc.col=rgb(0,.5,0,0.7), dist.col=rgb(0,0,0,0.3), sd.col=rgb(0,0,0,0.5), rule=1) {
  # set the calibration curve
  ccdir <- .validateDirectoryName(ccdir)
  if(ccdir == "")
    ccdir = paste(system.file("extdata", package="IntCal"), "/", sep="")

  # set calibration curve
  if(cc==1) calcurve <- read.table(paste(ccdir, cc1,  sep="")) else
    if(cc==2) calcurve <- read.table(paste(ccdir, cc2,  sep="")) else
      if(cc==3) calcurve <- read.table(paste(ccdir, cc3,  sep="")) else
        if(cc==4) calcurve <- read.table(paste(ccdir, cc4,  sep="")) else
            stop("I do not understand which calibration curve you mean, please check the manual", call.=FALSE)     

  # include postbomb curve if required
  if(cage < 0) {
    pb <- 0
    if(postbomb==FALSE)
      stop("\n  Negative 14C age, should I use a postbomb curve?\n", call.=FALSE)
    if(postbomb==1) pb <- pb1 else
       if(postbomb==2) pb <- pb2 else
         if(postbomb==3) pb <- pb3 else
           if(postbomb==4) pb <- pb4 else
             if(postbomb==5) pb <- pb5 else
               stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
	yrsteps <- min(pbsteps, yrsteps)
    if(length(pb) > 0) {
      pb <- read.table(system.file("extdata", pb, package=packageName()))
      pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
      pb.y <- approx(pb[,1], pb[,2], pb.x, rule=rule)$y
      pb.sd <- approx(pb[,1], pb[,3], pb.x, rule=rule)$y
      calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
    }
    cat("  postbomb date, interpolating to every", pbsteps, "yr.")
  }

  # check whether date lies partly or entirely beyond the calibration curve
  if(length(reservoir) == 2) { # assuming that first value is mean offset, second is error
    error <- sqrt(error^2 + reservoir[2]^2)
    reservoir <- reservoir[1]
  }
  border <- 0
  if(cage-reservoir-error < min(calcurve[,2]+calcurve[,3]))
    if(cage-reservoir+error > min(calcurve[,2]-calcurve[,3]))
      border <- 1 else border <- 2
  if(cage-reservoir+error > max(calcurve[,2]-calcurve[,3]))
    if(cage-reservoir-error < max(calcurve[,2]+calcurve[,3]))
      border <- 1 else border <- 2
  if(border == 1)
    cat("\nDate falls partly beyond calibration curve and will be truncated!")
  if(border == 2)
    stop("\nCannot calibrate dates beyond calibration curve!\n\n")

  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD) {
    theta <- 1950-calcurve[,1]
    ad <- max(which(theta > 0)) # one side of the border between AD and BC
    theta <- c(theta[1:(ad-1)], theta[ad]:theta[ad+2], theta[(ad+3):length(theta)])
   	mu <- approx(1950-calcurve[,1], calcurve[,2], theta, rule=rule)$y
    sigma <- approx(1950-calcurve[,1], calcurve[,3], theta, rule=rule)$y
    theta[theta <= 0] <- theta[theta <= 0] - 1
    calcurve <- cbind(theta, mu, sigma)
  } else 
      theta <- calcurve[,1]
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu
  f.cage <- exp(-(cage-reservoir)/8033)
  f.error <- f.cage - exp(-(cage-reservoir+error)/8033)

  # calibrate the date and report its highest posterior density (hpd) range
  if(length(xlab) == 0)
  xlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  calib <- .caldist(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD, rule=rule)
  hpd <- .hpd(calib, prob, hpdsteps, yrsteps)
  colnames(hpd) <- c("yrmin", "yrmax", "prob")
  dat <- list(calib=calib, hpd=hpd)
  if(storedat)
    dat <<- dat
  cat("\nmin\tmax\tprob\n")
  for(i in 1:nrow(hpd)) {
    for(j in 1:3)
	  cat(hpd[i,j], "\t")
    cat("\n")
  }
  cat("\n")

  # produce a graph of the calibrated distribution (default)
  if(graph) {
    ifelse(BCAD,
      xrange <- 1950+c((1+expand)*(min(calib[,1])-1950), (1-expand)*(max(calib[,1])-1950)),
      xrange <- c((1+expand)*max(calib[,1]), (1-expand)*min(calib[,1])))
    if(length(yrmin) > 0)
	  xrange[2] <- yrmin
    if(length(yrmax) > 0)
	  xrange[1] <- yrmax
    ifelse(BCAD,
      cc <- calcurve[max(which(theta >= min(xrange))):min(which(theta <= max(xrange))),],
      cc <- calcurve[min(which(theta >= min(xrange))):max(which(theta <= max(xrange))),])

  # first plot the calibrated distribution, and its hpd ranges
  par(mar=mar, mgp=mgp, bty=bty, xaxs=xaxs, xaxt="s", yaxs=yaxs, yaxt="n", new=FALSE)
  pol <- cbind(c(calib[,1], rev(calib[,1])), c(calib[,2]/max(calib[,2]), rep(0, length=nrow(calib))))
  plot(0, type="n", xlim=xrange, ylim=c(0,1/calheight), xlab="", ylab="")
  polygon(pol, col=dist.col, border=NA)
  for(i in 1:nrow(hpd)) {
    if(hpd[i,1]==hpd[i,2]) {
      probs <- calib[which(calib[,1]==hpd[i,1]),]
      lines(rep(probs[1], 2), c(0, probs[2]/max(calib[,2])), col=grey(.5))
    } else {
        probs <- calib[max(which(calib[,1]<=hpd[i,1])):max(which(calib[,1]<=hpd[i,2])),]
        pol <- cbind(c(probs[,1], rev(probs[,1])), c(probs[,2]/max(calib[,2]), rep(0, length=nrow(probs))))
        polygon(pol, col=sd.col, border=NA)
      }
  }
  lines(calib[,1], calib[,2]/max(calib[,2]))
  abline(h=0)

  # now draw the 14C distribution (normal distribution, on vertical axis)
  par(new=TRUE, yaxt="s", yaxs="r", xaxt="n")
  if(length(cc) == 3)
	cc <- cbind(cc[1], cc[2], cc[3])
  if(reservoir != 0)
    main <- substitute(cage-res %+-% er, list(cage=cage, er=error, res=reservoir)) else
      main <- substitute(cage %+-% er, list(cage=cage, er=error))
  if(length(title)>0)
	main <- title
  if(length(minC14) == 0)
    minC14 <- min(cc[,2]-qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir-qnorm(1-(1-prob)/2)*error)
  if(length(maxC14) == 0)
    maxC14 <- max(cc[,2]+qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir+qnorm(1-(1-prob)/2)*error)
  if(length(ylab) == 0)
    ylab <- expression(paste(""^14, "C BP"))
  plot(0, type="n", xlim=xrange, ylim=c(minC14, maxC14), xlab=xlab, ylab=ylab, main=main)
  if(length(calibt) > 0)
    times <- 5*times
  yage <- (cage-reservoir-times*error):(cage-reservoir+times*error) # must not be on F14C for plot
  if(length(calibt) < 2)
    xage <- dnorm(exp(-yage/8033), f.cage, f.error) else
      xage <- (calibt[2] + ((f.cage-exp(-yage/8033))^2) / (2*(f.error^2))) ^ -(calibt[1]+0.5)
  xage.plot <- xrange[1]-((xrange[1]-xrange[2])*calheight)*xage/max(xage)
  pol <- cbind(c(xage.plot, rep(xrange[1], length(xage))), c(yage, rev(yage)))
  polygon(pol, col=dist.col, border="black")

  # draw the highest posterior density (hpd) range of the 14C date
  xage[which(cumsum(xage)/sum(xage) > 1 - (1-prob)/2)] <- 0
  xage[which(cumsum(xage)/sum(xage) < (1-prob)/2)] <- 0
  xage <- xrange[1]-((xrange[1]-xrange[2])*calheight)*(xage/max(xage))
  pol <- cbind(c(xage, rep(xrange[1], length=length(xage))), c(yage, rev(yage)))
  polygon(pol, col=sd.col, border=FALSE)

  # plot the mid and error of the 14C date
  points(xrange[1]-.01*(xrange[1]-xrange[2]), cage, pch=19, col=date.col)
  lines(rep(xrange[1]-.01*(xrange[1]-xrange[2]), 2), c(cage-error, cage+error), lwd=2, col=date.col)

  # now draw the calibration curve
  pol <- cbind(c(theta, rev(theta)),
    c(calcurve[,2]-qnorm(1-(1-prob)/2)*calcurve[,3], rev(calcurve[,2]+qnorm(1-(1-prob)/2)*calcurve[,3])))
  polygon(pol, border=cc.col, col=cc.col)
  }
  invisible(dat)
}



#' @name mix.calibrationcurves
#' @title Build a custom-made, mixed calibration curve. 
#' @description If two curves need to be 'mixed' to calibrate, e.g. for dates of mixed terrestrial and marine carbon sources, then this function can be used. 
#' @details The proportional contribution of each of both calibration curves has to be set. 
#'
#' @param proportion Proportion of the first calibration curve required. e.g., change to \code{proportion=0.7} if \code{cc1} should contribute 70\% (and \code{cc2} 30\%) to the mixed curve.
#' @param cc1 The first calibration curve to be mixed. Defaults to the northern hemisphere terrestrial curve IntCal20.
#' @param cc2 The second calibration curve to be mixed. Defaults to the marine curve Marine20.
#' @param name Name of the new calibration curve.
#' @param dirname Directory where the file will be written. If using the default \code{dirname="."}, 
#' the new curve will be saved in current working directory. 
#' @param offset Any offset and error to be applied to \code{cc2} (default 0 +- 0).
#' @param sep Separator between fields (tab by default, "\\t") 
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A file containing the custom-made calibration curve, based on calibration curves \code{cc1} and \code{cc2}.
#' @examples
#'   mix.calibrationcurves(dirname=tempdir())
#' @export
mix.calibrationcurves <- function(proportion=.5, cc1="3Col_intcal20.14C", cc2="3Col_marine20.14C", name="mixed.14C", dirname=".", offset=c(0,0), sep="\t", rule=1) {
  ccloc <- normalizePath(system.file("extdata/", package='IntCal')) 
  dirname <- .validateDirectoryName(dirname)
  
  cc1 <- read.table(normalizePath(paste(ccloc, "/", cc1,  sep="")))
  cc2 <- read.table(normalizePath(paste(ccloc, "/", cc2,  sep="")))
  cc2.mu <- approx(cc2[,1], cc2[,2], cc1[,1], rule=rule)$y + offset[1] # interpolate cc2 to the calendar years of cc1
  cc2.error <- approx(cc2[,1], cc2[,3], cc1[,1], rule=rule)$y
  cc2.error <- sqrt(cc2.error^2 + offset[2]^2)
  mu <- proportion * cc1[,2] + (1-proportion) * cc2.mu
  error <- proportion * cc1[,3] + (1-proportion) * cc2.error
  write.table(cbind(cc1[,1], mu, error), paste(dirname, name,  sep="") , row.names=FALSE, col.names=FALSE, sep=sep)
}







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
student.t <- function(y=2450, error=50, t.a=3, t.b=4, cc=1, postbomb=NULL, cc1="IntCal20", cc2="Marine20", cc3="SHCal20", cc4="mixed", ccdir="",Cutoff=1e-5, times=8, rule=1)
{
  ccdir <-.validateDirectoryName(ccdir)
  # set the calibration curve
  if(ccdir=="")
    ccdir = paste(system.file("extdata", package="IntCal"), "/", sep="")

  if(cc == 0)
  {
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
  } else
  {
    
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
    else 
    {
      
      if(postbomb==1) bomb <- read.table(system.file("extdata","postbomb_NH1.14C", package=packageName()))[,1:3] else
        if(postbomb==2) bomb <- read.table(system.file("extdata","postbomb_NH2.14C", package=packageName()))[,1:3] else
          if(postbomb==3) bomb <- read.table(system.file("extdata","postbomb_NH3.14C", package=packageName()))[,1:3] else
            if(postbomb==4) bomb <- read.table(system.file("extdata","postbomb_SH1-2.14C", package=packageName()))[,1:3] else
              if(postbomb==5) bomb <- read.table(system.file("extdata","postbomb_SH3.14C", package=packageName()))[,1:3] else
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


# find the 14C age of a single cal BP year
#' @name calBP.14C
#' @title Find the 14C age and error belonging to a cal BP age.
#' @description From a given calendar ages, the calibration curve (default cc=1) is interpolated and the corresponding 14C age and error is returned.
#' @details Interpolation is used, and values outside the calibration curve are given as NA. For negative cal BP ages, a postbomb curve will have to be provided. 
#' @param yr The cal BP year.
#' @param cc calibration curve for C14 dates (1, 2 or 3).
#' @param cc1 For northern hemisphere terrestrial C14 dates.
#' @param cc2 For marine C14 dates.
#' @param cc3 For southern hemisphere C14 dates.
#' @param cc4 A custom calibration curve
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are allocated. By default \code{ccdir=""}. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param postbomb Which postbomb curve to use for negative 14C dates
#' @param pb1 For Northern hemisphere region 1 postbomb C-14 dates.
#' @param pb2 For Northern hemisphere region 2 postbomb C-14 dates.
#' @param pb3 For Northern hemisphere region 3 postbomb C-14 dates.
#' @param pb4 For Southern hemisphere regions 1-2 postbomb C-14 dates.
#' @param pb5 For Southern hemisphere region 3 postbomb C-14 dates.
#' @param rule How should R's approx function deal with extrapolation. If \code{rule=1}, the default, then NAs are returned for such points and if it is 2, the value at the closest data extreme is used.
#' @author Maarten Blaauw
#' @examples 
#' calBP.14C(100) 
#' 
#' @export
calBP.14C <- function(yr, cc=1, cc1="3Col_intcal20.14C", cc2="3Col_marine20.14C", cc3="3Col_shcal20.14C", cc4="mixed.14C", ccdir="", postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C", pb5="postbomb_SH3.14C", rule=1) {
  # set the calibration curve
  ccdir <- .validateDirectoryName(ccdir)
  if(ccdir == "")
    ccdir = paste(system.file("extdata", package="IntCal"), "/", sep="")

  # set calibration curve
  if(cc==1) calcurve <- read.table(paste(ccdir, cc1,  sep="")) else
    if(cc==2) calcurve <- read.table(paste(ccdir, cc2,  sep="")) else
      if(cc==3) calcurve <- read.table(paste(ccdir, cc3,  sep="")) else
        if(cc==4) calcurve <- read.table(paste(ccdir, cc4,  sep="")) else
           stop("I do not understand which calibration curve you mean, please check the manual", call.=FALSE)     

  # use postbomb curve if required
  if(yr < 0) {
    pb <- 0
    if(postbomb==FALSE)
      stop("\n  Negative 14C age, should I use a postbomb curve?\n", call.=FALSE)
    if(postbomb==1) pb <- pb1 else
       if(postbomb==2) pb <- pb2 else
         if(postbomb==3) pb <- pb3 else
           if(postbomb==4) pb <- pb4 else
             if(postbomb==5) pb <- pb5 else
               stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
    if(length(pb) > 0) 
      calcurve <- read.table(system.file("extdata", pb, package=packageName()))	
  }
  mu <- approx(calcurve[,1], calcurve[,2], yr, rule=rule)$y
  er <- approx(calcurve[,1], calcurve[,3], yr, rule=rule)$y
  return(c(mu, er))
}
 

# find the calibrated distributions of 14C dates
.caldist <- function(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD, normalise=FALSE, rule=rule)
  {
    if(f.cage > 1)
      {
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


# find the highest posterior density (hpd) of the calibrated distribution
.hpd <- function(dat, prob, hpdsteps, yrsteps, rule=1)
  {
    # interpolate and rank the ages according to their calibrated distribution probabilities
    dat <- approx(dat[,1], dat[,2], seq(min(dat[,1]), max(dat[,1]), by=yrsteps), rule=rule)
    o <- order(dat$y, decreasing=TRUE)
    dat <- cbind(dat$x[o], dat$y[o]/sum(dat$y))

    # only retain those ages with cumulative normalised probabilities within required percentage
    dat <- dat[which(cumsum(dat[,2]) <= prob),]
    dat <- dat[order(dat[,1]),]

    # identify any individual ranges within the hpd range and calculate their probability
    dif <- which(diff(dat[,1]) > hpdsteps)
    if(length(dif)==0)
      hpds <- cbind(min(dat[,1]), max(dat[,1]), 100*prob) else
        {
          dif <- c(dat[1,1], sort(c(dat[dif,1], dat[dif+1,1])), dat[nrow(dat),1])
          dif <- matrix(dif, ncol=2, byrow=TRUE)
          #probs <- c()
          probs <- numeric(nrow(dif))
          for(i in 1:nrow(dif))
            probs[i] <- round(100*sum(dat[which(dat[,1]==dif[i,1]):which(dat[,1]==dif[i,2]),2]), 1)
          hpds <- cbind(dif, probs)
        }
    hpds
  }


