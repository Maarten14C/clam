
#' clam
#' 
#' clam is R code for classical (non-Bayesian) age-depth modelling.
#' 
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @importFrom grDevices dev.off grey rgb pdf png
#' @importFrom graphics abline image layout legend lines par plot points polygon rect
#' @importFrom stats approx density dnorm lm loess pnorm predict qnorm quantile rnorm runif smooth.spline spline weighted.mean
#' @importFrom utils read.csv read.table write.table packageName
#' @name clam
NULL  

# done: 

# do: 

#' @name clam
#' @title The main age-depth modelling function
#' @description Produce age-depth models for cores with dated depths.
#' @details 
#' Cores containing several 14C and/or other dates can be processed semi-automatically in order to obtain age-depth models. 
#' In the process, any 14C dates are calibrated, and age-depth curves are repeatedly drawn through point estimates sampled from the dates.
#'  Age-depth models can be based on linear interpolation, linear/polynomial regression, or cubic, smooth or locally weighted splines. 
#'  For each date, the probability of a calendar year being sampled is proportionate to its calibrated probability (see Blaauw, 2010). 
#'  Uncertainty ranges as well as a 'best' age-model are calculated.
#'
#' Additional cores should be put in a comma-separated file in a sub-folder of the directory where the cores are stored. 
#' By default this parent folder is called \code{coredir="clam_runs"} (if no folder called \code{"Cores"} already exists). If your core is called MyCore1, save MyCore1.csv as \code{clam_runs/MyCore1/MyCore1.csv}.
#' Ensure that the names of the core's folder and filename's root (the part before .csv) match, e.g., using exactly similar upper- and lower case letters.
#' 
#' Avoid the use of spaces or non-standard (non-ASCII) characters within the file or in folder or file names. 
#' The plain text file should consist of 6 or 7 columns (also called fields), containing in the following exact order (see the example below):
#'  \enumerate{
#'   \item Identification labels (e.g. 14C lab codes)
#'   \item 14C ages for 14C-dated depths; leave empty for non-14C dated depths
#'   \item cal BP ages (for any non-14C dates such as the core surface; leave empty for levels with 14C dates)
#'   \item errors (reported 1 standard deviation errors. This column should never be left empty. Errors should always be larger than 0)
#'   \item age offsets if known (otherwise leave empty)
#'   \item depths (depths in the sequence were the dated samples were taken, default unit depth="cm"; this column should never be left empty)
#'   \item thicknesses of the sampled slices (optional column; leave empty for default of 1)
#'  }
#' Add a final empty line to your core's .csv file by pressing 'Enter' after the file's last value. 
#' 
#' These files can be made in spreadsheet software such as MS-Excel, but it is always a good idea to check the file's formatting in a plain-text editor such as WordPad. Remove any lines which contain only commas, and it is also recommended to remove quotes ()\code{\" or \'}) in the headers or elsewhere. 
#'
#' Age-models for the core can then be produced by typing, e.g., \code{clam("MyCore1")}.
#'
#' By default the northern hemisphere terrestrial calibration curve is used (\code{cc=1}, \code{cc1="IntCal20.14C"}). 
#' To use alternative curves, change \code{cc} to \code{cc=2 (cc2="Marine20.14C")}, \code{cc=3 (cc3="SHCal20.14C")}, \code{cc=4 (cc4="mixed.14C")}. 
#' You can also provide custom-built calibration curves, indicating its location using \code{ccdir}.
#'
#' The provided example (default \code{core="Example"}) is core Quilichao-1 which was sampled from a Colombian lake (Berrio et al., 2002). 
#' This core was chosen because it was dated at a rather high resolution, and appears to contain a hiatus (e.g., try \code{hiatus=450}
#' for a hiatus at 450 cm depth).
#'
#' Each clam run will produce a range of files within the core's folder. One, ending with \code{"_calibrated.txt"} contains the calibrated 
#' age ranges of the 14C and other dates. The others will be named according to the core's name followed by the model type, 
#' and contain the age estimates for all depths (files ending with \code{"_ages.txt"}), settings (files ending with \code{"_settings.txt"}) 
#' and graphs (files ending with \code{".pdf"} and \code{".png"}). 
#' The file containing the age estimates has 5 columns; first the depths, then the minima and maxima of the confidence intervals, 
#' then a "best" estimate, and finally the reconstructed accumulation rates. The reported values are rounded to 0 decimals by default
#' (\code{decimals=0}). Accumulation rates are in yr/cm ("deposition time") by default (\code{cmyr=FALSE}), but can be reported in cm/yr (\code{cmyr=TRUE}).
#'
#' see accompanying webpage \url{http://www.qub.ac.uk/chrono/blaauw/clam.html} and Blaauw 2010 (Quaternary Geochronology 5: 512-518).
#'
#' @param core Name of the core, given using quotes. Defaults to the core provided with clam, \code{core="Example"}.
#' @param type The type of age-depth model. Five different types are provided:
#'  \enumerate{
#'   \item linear interpolation between neighbouring levels (1, "int", "inter" or "interp") 
#'   \item linear or higher polynomial regression (2, "reg", "regr", "poly" or "polyn", default linear) 
#'   \item cubic spline (3, "spl" or "spline") 
#'   \item smooth spline (4, "sm" or "smooth", default smoothing 0.3) 
#'   \item locally weighted spline (5, "loess" or "lowess", default smoothing 0.75, cannot extrapolate)
#'  }
#' @param smooth Degree of smoothing. Gives polynomial degree for model type 2. Not relevant for \code{type=1 or type=3}.
#'  \itemize{
#'   \item for type=2: \code{smooth=1} (linear), \code{smooth=2} second-order polynomial, \code{smooth=3} for third-order polynomial, etc. 
#'   \item for type=4: \code{smooth=0.3} 
#'   \item for type=5: \code{smooth=0.75} 
#'  }
#' @param prob Confidence intervals (between 0 and 1), default \code{prob=0.95} or 95\%.
#' @param its Amount of age-model iterations; defaults to \code{its=1000}.
#' @param coredir The directory where core runs are stored (each core in its own directory named after the core's name). 
#' @param ask By default, and as per R rules, clam will ask if it is OK to make or write to a directory.
#' Defaults to \code{coredir="clam_runs"}, or to \code{coredir="Cores"} if this folder exists where R is working. 
#' @param wghts Weights can be applied to dated depths as follows:
#'  \itemize{
#'   \item 0 no weighting 
#'   \item 1 weighted to calibrated probabilities of sampled calendar years (default, \code{wghts=1}). 
#'   \item 2 weighted to (inverse squared) errors of the dates.
#'  }
#' @param cc calibration curve for C14 dates (1, 2 or 3).
#' @param cc1 For terrestrial, northern hemisphere C14 dates.
#' @param cc2 For marine C14 dates.
#' @param cc3 For southern hemisphere C14 dates.
#' @param cc4 For mixed terrestrial/marine C14 dates.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}. See \url{http://calib.org/CALIBomb/}.
#' @param pb1 For Northern hemisphere region 1 postbomb C14 dates.
#' @param pb2 For Northern hemisphere region 2 postbomb C14 dates.
#' @param pb3 For Northern hemisphere region 3 postbomb C14 dates.
#' @param pb4 For Southern hemisphere regions 1-2 postbomb C14 dates.
#' @param pb5 For Southern hemisphere region 3 postbomb C14 dates.
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""}. 
#' For example, use \code{ccdir="."} to choose current working directory, or \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param outliers The number of any dates to be considered outlying, e.g. \code{c(5,6)} for the fifth and sixth dated depth counting from the top of a core.
#' @param ignore The number of any dates that should be ignored, e.g., \code{c(5,6)} for the fifth and sixth date counting from the top of a core.
#' @param youngest The age beyond which dates should be truncated (e.g., \code{youngest=-60} if the core was sampled in -60 cal BP or AD 2010).
#' @param extradates Depths of any additional dates with their files of ages and probabilities.
#' @param slump Upper and lower depths of sections of abrupt accumulation that should be excised, e.g., \code{c(600, 550, 120, 100)} for two sections of 600-550 and 120-100 cm depth.
#' @param est Which point estimate to use as 'best' age. It is highly recommended to not only use these 'best' point estimates, as chronological uncertainties are often considerable and should not be ignored.  
#'  \enumerate{
#'   \item averages of age-depth model derived ages (default, \code{est=1})
#'   \item midpoints of age-depth model derived age estimates
#'   \item midpoints of calibrated ranges
#'   \item weighted means of calibrated ranges 
#'   \item medians of calibrated distributions 
#'   \item maximum densities of calibrated distributions 
#'   \item midpoints of entire calibrated distributions (including years outside the calibrated ranges)
#'  }
#' @param calibt Calibration based on the student-t distribution. By default, the Gaussian distribution is used (\code{calibt=FALSE}). To use the student-t distribution, provide two parameters such as \code{calibt=c(3,4)}.
#' @param mixed.effect Set to \code{TRUE} to activate mixed-effect modelling.
#' @param dmin Minimum depth of age-depth model (e.g., extrapolate).
#' @param dmax Maximum depth of age-depth model (e.g., extrapolate).
#' @param every Resolution at which (ages for) depths are calculated.
#' @param yrmin Minimum of calendar axis of age-depth plot (calculate automatically by default).
#' @param yrmax Maximum of calendar axis of age-depth plot (calculated automatically by default).
#' @param yrsteps Temporal resolution at which calibrated ages are calculated (in calendar years).
#' @param pbsteps Temporal resolution at which postbomb C14 ages are calibrated (in calendar years).
#' @param hpdsteps Temporal resolution at which highest posterior density ranges are calibrated (in calendar years).
#' @param BCAD Use BC/AD or cal BP scale.
#' @param decimals Amount of decimals for rounding.
#' @param cmyr Accumulation rates can be provided as yr/cm (default, \code{cmyr=TRUE}, more accurately named deposition times) or cm/yr (\code{cmyr=FALSE}).
#' @param ageofdepth Calculate age estimates of a specific depth.
#' @param depth Depth units.
#' @param depthseq Sequence of depths for which age estimates are to be calculated (default: from \code{dmin} to \code{dmax} with steps of size every)
#' @param depths.file Use a file with depths for depthseq.
#' @param thickness Thickness of the dated samples.
#' @param hiatus Depths of any hiatuses, e.g., \code{c(500, 300)}. Each sub-section must have at least 2 dates (\code{4} for smoothing spline; does not work with loess as it cannot extrapolate).
#' @param remove.reverse Proportion of age-models with reversals that can be removed before prompting a warning.
#' Set at \code{FALSE} to avoid removing models with reversals.
#' @param times Half-range of calibration curve used to calibrate dates (multiplication factor for the dates' errors).
#' @param sep Separator between the fields of the plain text file containing the dating information.
#' @param ext Extension of the file containing the dating information.
#' @param runname Text to add to the core name for specific runs, e.g., "MyCore_Test1"
#' @param storedat Store the dates and age-model within R after a \code{clam} run. Defaults to \code{storedat=TRUE}.
#' @param threshold Below which value should probabilities be excluded from calculations.
#' @param proxies Set to \code{TRUE} to plot proxies against age after the run.
#' @param revaxes Set to \code{TRUE} to plot ages on the vertical axis and depth on the horizontal axis.
#' @param revd Plot depth axis in reverse.
#' @param revyr Plot age axis in reverse.
#' @param calhght Heights of the calibrated distributions in the age-depth plot.
#' @param maxhght Maximum height of age probability distributions.
#' @param mirror Plot the age distributions in "mirror" style (above and below depth).
#' @param plotrange Plot the confidence ranges of the age-model.
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param bty Type of box to be drawn around plots. Draw a box around the graph (\code{"n"} for none, 
#' and \code{"l"}, \code{"7"}, \code{"c"}, \code{"u"}, "]" or \code{"o"} for correspondingly shaped boxes).
#' @param plotpdf Produce a pdf file of the age-depth plot.
#' @param plotpng Produce a png file of the age-depth plot.
#' @param greyscale Produce a grey-scale representation of all age-models (number gives resolution, e.g., 500 bins; will cancel plotting of the confidence intervals).
#' @param yrlab Label of the calendar axis. Defaults to either cal BP or BC/AD. Alternative names can be provided.
#' @param dlab Label of the depth axis. Defaults to \code{dlab="Depth (cm)"} (assuming \code{depth="cm"}), but alternative names can be provided.
#' @param calcol Colour of the calibrated distributions in the age-depth plot.
#' @param C14col Colour of the calibrated ranges of the dates.
#' @param outcol Colour of outlying dates.
#' @param outlsize Size of symbols outlying dates.
#' @param bestcol Colour of the "best" age-depth model (based on chosen value for est).
#' @param rangecol Colour of plotted confidence ranges.
#' @param slumpcol Colour of slump.
#' @param plotname Print the core name on the graph.
#' @param ash Plot all distributions at the same height.
#' 
#' @author Maarten Blaauw
#' @return Age model construction together with a text output and files saved to a folder in the \code{coredir/core} directory.
#' @examples 
#'  clam(, coredir=tempdir()) # Create the example in Cores/Example folder
#'  clam(, coredir=tempdir(), extradates=470) 
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/clam.html}
#' \link{calibrate}
#' \link{mix.calibrationcurves}
#' \link{pMC.age}
#' \link{age.pMC}
#' \link{student.t}
#' \link{deptime.depth}
#' \link{deptime.age}
#' \link{plot_proxies}
#' 
#' @references
#' Berrio, J.C., Hooghiemstra, H., Marchant, R., Rangel, O., 2002. Late-glacial and Holocene history of the dry forest area 
#' in the south Colombian Cauca Valley. _Journal of Quaternary Science_ 17, 667-682
#'
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
clam <- function(core="Example", type=1, smooth=NULL, prob=0.95, its=1000, coredir=NULL, ask=TRUE, wghts=1, cc=1, cc1="IntCal20.14C", cc2="Marine20.14C", cc3="SHCal20.14C", cc4="mixed.14C",  postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C",pb5="postbomb_SH3.14C", ccdir="", outliers=NULL, ignore=NULL, youngest=NULL, extradates=NULL, slump=NULL, est=1, calibt=FALSE, mixed.effect=FALSE, dmin=NULL, dmax=NULL, every=1, yrmin=NULL, yrmax=NULL, yrsteps=1, pbsteps=0.01, hpdsteps=1, BCAD=FALSE, decimals=0, cmyr=FALSE, ageofdepth=NULL, depth="cm", depthseq=NULL, depths.file=FALSE, thickness=1, hiatus=NULL, remove.reverse=0.5, times=5, sep=",", ext=".csv", runname=NULL, storedat=TRUE, threshold=1e-6, proxies=FALSE, revaxes=FALSE, revd=TRUE, revyr=TRUE, calhght=0.3, maxhght=0.01, mirror=TRUE, plotrange=TRUE, bty="l", mar=c(3.5,3,2,1), mgp=c(2,1,0), plotpdf=TRUE, plotpng=TRUE, greyscale=NULL, yrlab=NULL, dlab=NULL, calcol=rgb(0,0.5,0.5,0.5), C14col=rgb(0,0,1,0.5), outcol="red", outlsize=1, bestcol="black", rangecol=rgb(0,0,0,0.3), slumpcol=grey(0.75), plotname=TRUE, ash=FALSE) {
  # If coredir is left empty, check for a folder named Cores in the current working directory, and if this doesn't exist, for a folder called clam_runs (make this folder if it doesn't exist yet).
  # Check if we have write access. If not, tell the user to provide a different, writeable location for coredir. 

  coredir <- assign_coredir(coredir, core, ask)	  
  if(core == "Example") {
    dir.create(paste(coredir, "Example/", sep=""), showWarnings = FALSE, recursive = TRUE)
    fileCopy <- system.file("extdata/Example/", package="clam")
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite=FALSE)
  } 
 
  # set the calibration curve
  if(ccdir == "")
    ccdir <- system.file("extdata", package=packageName()) 

  # warn and stop if unexpected settings are provided
  if(type > 5 || type < 1 || prob < 0 || prob > 1 || its < 100 || wghts < 0 || wghts > 1 || est < 1 || est > 7 || yrsteps <= 0 || hpdsteps <= 0 || every <= 0 || decimals < 0 || cmyr < 0 || cmyr > 1 || thickness < 0 || times < 1 || calhght < 0 || (type==5 && length(hiatus)>0))
    stop("\n Warning, clam cannot run with these settings! Please check the manual.\n\n", call.=FALSE)

  csvFile <- paste(coredir, core, "/", core, ext, sep="")
  if(!file.exists(csvFile))
    stop(paste("\nInput data file", csvFile, "not found.", sep=" "), call.=FALSE)
  dets <- read.csv(csvFile, sep=sep)
  d <- dets[,6]
  if(min(diff(d)) < 0)
    cat("\n Warning, depths not in ascending order (top ones should come first).\n\n")

  # avoid Windows/Mac habit of silently adding .txt extension to plain text files
  Gates <- list.files(paste(coredir, core, sep=""), pattern=".csv.txt")
  if(length(Gates) > 0) {
    cat("\nRemoving unnecessary .txt extension from .csv file", Gates[1], "\n")
    file.rename(paste(coredir, core, "/", core, ".csv.txt", sep=""),
      paste(coredir, core, "/", core, ".csv", sep=""))
  }

  # set the calibration curve
  ccdir <- .validateDirectoryName(ccdir)
  if(ccdir == "") # so, if no alternative folder provided, use clam's calibration curves
    ccdir = paste(system.file("extdata", package=packageName()), "/", sep="")
    
  if(cc==1) calcurve <- read.table(paste(ccdir, cc1,  sep="")) else
    if(cc==2) calcurve <- read.table(paste(ccdir, cc2,  sep="")) else
      if(cc==3) calcurve <- read.table(paste(ccdir, cc3,  sep="")) else
        if(cc==4) calcurve <- read.table(paste(ccdir, cc4,  sep="")) else
          stop("I do not understand which calibration curve you mean, check the manual", call.=FALSE)
  if(cc==1) ccname <- cc1 else
    if(cc==2) ccname <- cc2 else
      if(cc==3) ccname <- cc3 else
        if(cc==4) ccname <- cc4 

  # negative C14 ages need a postbomb curve
  pbnames <- c(pb1, pb2, pb3, pb4, pb5)
  cdat <- dets[,2]
  if(length(cdat[!is.na(cdat)]) > 0)
    if(min(cdat[!is.na(cdat)]) < 0)
      if(postbomb==FALSE)
        cat("Warning, negative 14C ages, should I use a postbomb curve?") else {
          if(postbomb>5)
            stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
          yrsteps <- min(pbsteps, yrsteps)
          pb <- read.table(system.file("extdata", pbnames[postbomb], package=packageName()))
          pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
          pb.y <- approx(pb[,1], pb[,2], pb.x)$y
          pb.sd <- approx(pb[,1], pb[,3], pb.x)$y
          calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
        }

  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD) {
    theta <- 1950-calcurve[,1]
    border <- max(which(theta > 0))
    theta <- c(theta[1:(border-1)], theta[border]:theta[border+2], theta[(border+3):length(theta)])
    mu <- approx(1950-calcurve[,1], calcurve[,2], theta)$y
    sigma <- approx(1950-calcurve[,1], calcurve[,3], theta)$y
    theta[theta <=0] <- theta[theta <= 0]-1
    calcurve <- cbind(theta, mu, sigma)
  } else theta <- calcurve[,1]
  if(length(yrlab) == 0)
    yrlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu

  # prepare for slumps and hiatuses
  #cat(paste("Core name Test:", core))
  if(length(greyscale) > 0) storedat <- TRUE
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("\n Warning, slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
    if(length(dmax)==0)
      dmax <- max(dets[,6])
    if(length(extradates) > 0)
      dmax <- max(dmax, extradates)
    for(i in 1:nrow(slump)) {
      d[d > min(slump[i,])] <- d[d > min(slump[i,])] - (max(slump[i,]) - min(slump[i,]))
      dmax <- dmax - (max(slump[i,])-min(slump[i,]))
    }
  if(length(hiatus) > 0)
    for(i in 1:nrow(slump)) {
       below.slump <- which(hiatus > max(slump[i,]))
       above.slump <- which(hiatus < min(slump[i,]))
       hiatus[below.slump] <- hiatus[below.slump] - (max(slump[i,])-min(slump[i,]))
       hiatus <- hiatus[c(above.slump, below.slump)]
    }
  }

  # read in the data
  dat <- .read.clam(core, coredir, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb)
  cat("\n Calibrating dates... ")

  # calculate the depths to be used, based on the ranges and resolution
  if(length(dmin)==0)
    dmin <- floor(min(dat$depth))
  if(length(dmax)==0)
    dmax <- ceiling(max(dat$depth))
  if(depths.file)
    if(file.exists(dd <- paste(coredir, core, "/", core, "_depths.txt", sep=""))) {
      if(length(depthseq) == 0)
        depthseq <- seq(dmin, dmax, by=every)
      depthseq <- sort(unique(c(depthseq, suppressWarnings(read.table(dd))[,1])))
      dmin <- min(depthseq)#, read.table(dd)[,1])
      dmax <- max(depthseq)#, read.table(dd)[,1])
    } else
      stop(paste("\nCannot find file ", dat$core, "_depths.txt!\n", sep=""), call.=FALSE)
  if(length(depthseq) == 0)
    depthseq <- seq(dmin, dmax, by=every)
  if(proxies) {
    storedat <- TRUE
    if(file.exists(dd <- paste(coredir, core, "/", core, "_proxies.csv", sep="")))
      dat$proxies <- suppressWarnings(read.csv(dd, sep=sep)) else
        stop(paste("\nCannot find file ", dat$core, " _proxies.csv!\n", sep=""), call.=FALSE)
    dmin <- min(depthseq, dat$proxies[,1])
    dmax <- max(depthseq, dat$proxies[,1])
    depthseq <- sort(unique(c(depthseq, dat$proxies[,1])))
  }

  if(length(ageofdepth) > 0)
    depthseq <- sort(unique(c(ageofdepth, depthseq)))

  # decide which models and point estimates should be used
  if(any(type==c(1, "int", "inter", "interp"))) type <- 1 else
    if(any(type==c(2, "reg", "regr", "poly", "polyn"))) type <- 2 else
      if(any(type==c(3, "spline", "spl"))) type <- 3 else
        if(any(type==c(4, "smooth", "sm"))) type <- 4 else
          if(any(type==c(5, "loess", "lowess"))) type <- 5
  best <- cbind(dat$mid1, dat$mid1, dat$mid1, dat$wmn, dat$med, dat$mode, dat$mid2)
  Est <- best[,est]  

  # remove outliers from the age-depth modelling
  if(length(outliers) > 0) {
    depths <- dat$depth[-outliers]
    errors <- dat$error[-outliers]
    calibs <- dat$calib[-outliers]
    Est <- Est[-outliers]
  } else {
      depths <- dat$depth
      errors <- dat$error
      calibs <- dat$calib
    }

  # age-depth modelling with curves through sampled age estimates
  # in sections if one or more hiatuses are present
  if(length(hiatus) > 0) {
    allrange <- c(0,0,0,0)
    hiatusseq <- sort(c(range(depthseq), hiatus))
    for(i in 2:length(hiatusseq)) {
      cat(paste("\n section ", i-1, ",", sep=""))
      section <- depthseq[min(which(depthseq >= hiatusseq[i-1])) : max(which(depthseq <= hiatusseq[i]))]
      if(i>2) section <- section[-1]
      sel <- min(which(depths >= min(section))):max(which(depths <= max(section)))
      if(mixed.effect)
        if(length(outliers) > 0)
          smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) else
            smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt) else
              smp <- .smpl(its, depths[sel], calibs[sel], Est[sel])
      calrange <- .model.clam(type, smooth, its, wghts, depths[sel], errors[sel], section, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
      allrange <- rbind(allrange, calrange)
    }
  calrange <- allrange[2:nrow(allrange),]
  } else {
    if(mixed.effect)
	  if(length(outliers) > 0)
	    smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt) else
	      smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) else
	        smp <- .smpl(its, depths, calibs, Est)
     calrange <- .model.clam(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
    }
  dat$model <- approx(calrange[,1], (calrange[,2]+calrange[,3])/2, dat$depth)$y

  if(est==2) calrange[,4] <- (calrange[,2]+calrange[,3])/2
  if(!BCAD && any(diff(calrange[,4]) < 0) || BCAD && any(diff(calrange[,4]) > 0))
    reversal <- TRUE else reversal <- FALSE
  gfit <- round(.gfit(theta, f.mu, f.sigma, dat, calrange, outliers), 2)

  # re-correct the depths if slumps were applied
  if(length(slump) > 0) {
    dat <- .read.clam(core, coredir, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump=NULL, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb) # read in the original dates again
    calrange <- calrange[which(calrange[,1] <= dmax),]
    d <- calrange[,1]
    for(i in 1:nrow(slump)) {
      d[d > min(slump[i,])] <- d[d > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
      dmax <- dmax + (max(slump[i,]) - min(slump[i,]))
      calrange[,1] <- d
      hiatus[hiatus > min(slump[i,])] <- hiatus[hiatus > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
    }
  }

  # produce the age-depth plot, and a pdf copy if desired
  if(length(yrmin)==0)
    yrmin <- min(dat$mid1, calrange[,2])
  if(length(yrmax)==0)
    yrmax <- max(dat$mid1, calrange[,3])
  if(length(ageofdepth > 0))
    layout(matrix(c(1,2,1,3), nrow=2), heights=c(.7,.3))
  .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, yrlab, dlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) get('chron') else NULL, C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, core, bty, mar, mgp, ash)

  # write files providing calibrated dates, age-model and settings
  colnames(calrange) <- c("Depth", paste("min.", 100*prob, "%range", sep=""), paste("max.", 100*prob, "%range", sep=""), "point")
  .write.clam(dat, coredir, runname, calrange, core, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, cmyr, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, dlab, yrlab, plotrange, greyscale, if(length(greyscale)>0) get('chron') else NULL, C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
  closeAllConnections()

  if(storedat) {
    calrange <<- calrange
    dat <<- dat
    smp <<- smp
  }

  # plot the age distribution of a provided depth
  if(length(ageofdepth) > 0) {
    if(revaxes)
      abline(v=ageofdepth, lty=2) else
        abline(h=ageofdepth, lty=2)
    xlim <- range(.ageofdepth)
    if(!BCAD) xlim <- xlim[2:1]
    .ageofdepth <- get(".ageofdepth")
    hst <- density(.ageofdepth, n=max(1, max(xlim)-min(xlim)))
    yr <- seq(min(xlim), max(xlim), by=yrsteps)
    hst <- cbind(c(yr, max(xlim), min(xlim)), c(approx(hst$x, hst$y, yr)$y, 0, 0))
    plot(hst, type="n", main="", xlim=xlim, xlab=yrlab, ylab="")
    polygon(hst, col="grey")
    legend("topleft", paste(ageofdepth, depth), bty="n")
    layout(matrix(1))
    rng <- round(calrange[max(which(calrange[,1] <= ageofdepth)),])
    cat("\n  Age range of ", ageofdepth, " ", depth, ": ", rng[3], " to ", rng[2], ifelse(BCAD, " cal BC/AD", " cal BP"), " (", rng[3]-rng[2], " yr, ", prob, " % range)",  sep="")
  }

  # report the confidence ranges, the goodness-of-fit, and whether any age-reversals occurred
  rng <- round(calrange[,3]-calrange[,2])
  cat("\n  ", core, "'s ", 100*prob, "% confidence ranges span from ", min(rng), " to ", max(rng), " yr (average ", round(mean(rng)), " yr)", sep="")
  cat("\n  Fit (-log, lower is better):", gfit, "\n")
  if(reversal) 
    cat("  Age reversals occurred. Try other model?\n")
}



# calculate the age-depth model and its uncertainty
.model.clam <- function(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
  {
    # warn for extrapolation, refuse to do so for loess
    if(min(depthseq) < min(dat$depth) || max(depthseq) > max(dat$depth))
      if(type==5)
        stop(" cannot extrapolate using loess! Change settings.\n ", call.=FALSE) else
          cat(" extrapolating beyond dated levels, dangerous!\n ")

    # choose model: interpolation, (polynomial) regression, spline, smooth spline or loess
    chron <- array(0, dim=c(length(depthseq), its))
    if(type==1) chron <- .interp(depthseq, depths, its, chron, smp) else
      if(type==2) chron <- .poly(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
        if(type==3) chron <- .spline(depthseq, smooth, depths, its, chron, smp) else
          if(type==4) chron <- .smooth(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
            if(type==5) chron <- .loess(depthseq, smooth, wghts, errors, depths, its, chron, smp)

    # test against age reversals
    # warp <- c()
    warp <- NULL
    if(remove.reverse!=FALSE)
      for(i in 1:ncol(chron))
        if(!BCAD && min(diff(chron[,i])) <= 0 || BCAD && max(diff(chron[,i])) >= 0)
          warp <- c(warp, i)
    if(length(warp) > 0)
      if(length(warp) > remove.reverse*its)
        cat("\n\n !!! Too many models with age reversals!!!\n") else
          {
            cat("\n Removing", length(warp), "models with age reversals,", its-length(warp), "models left...")
            chron <- chron[,-warp]
            smp <- smp[,-warp,]
          }

    if(length(ageofdepth) > 0)
      if(ageofdepth %in% depthseq)
        .assign_to_global(".ageofdepth", chron[which(depthseq==ageofdepth),]) 

    if(storedat)
      {
        chron <<- chron
        smp <<- smp
      }

    # find uncertainty ranges of calendar age for each depth of the core
    calrange <- array(0, dim=c(nrow(chron), 2))
    # mn <- c()
    mn <- numeric(nrow(chron))
    for(i in 1:nrow(chron))
      {
        x <- chron[i,2:ncol(chron)]
        qp <- (1-prob)/2
        calrange[i,] <- quantile(x, c(qp, 1-qp))
        if(est==1) mn[i] <- mean(x)
      }
    if(est==1)
      cbind(depthseq, cbind(calrange, mn)) else
        cbind(depthseq, cbind(calrange, chron[,1]))
  }


