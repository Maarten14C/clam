
# read the data and perform first calculations incl. calibrations
.read.clam <- function(name, namedir,ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb)
  {
    coredir=paste(namedir, name, "/", sep="")
	if(!file.exists(paste(namedir, name, sep="")))
      stop(paste("\n\n Warning, cannot find a folder within", namedir," named ", name, ". Have you saved it in the right place and with the right name? Please check the manual\n\n", sep=""), call.=FALSE)
    if(!file.exists(paste(coredir, name, ext, sep="")))
      stop(paste(" \n\n Warning, cannot find file ", name, ".csv in folder",namedir, name, ". Have you saved it in the right place and named it correctly? Please check the manual\n\n", sep=""), call.=FALSE)
    dets <- suppressWarnings(read.table(paste(coredir, name, ext, sep=""), comment.char="", header=TRUE, sep=sep, na.strings = c("#N/A!", "NA", "@NA")))

    # read the file with the dating information
    dat <- list(coredir=coredir, name=name, calib=list(), ignore=NULL, ID=character(nrow(dets)), cage=numeric(nrow(dets)), 
    error=numeric(nrow(dets)), f.cage=numeric(nrow(dets)), f.error=numeric(nrow(dets)), outside=NULL, cal=NULL, res=NULL, depth=NULL, thick=NULL, BCAD=NULL,
    hpd=list(), mid1=numeric(nrow(dets)), mid2=numeric(nrow(dets)), wmn=numeric(nrow(dets)), med=numeric(nrow(dets)), mode=numeric(nrow(dets))) 

    # ignore dates if required, add thickness column if it was left out
    if(length(ignore) > 0)
      {
        dat$ignore <- as.character(dets[ignore,1])
        dets <- dets[-ignore,]
      }
    if(ncol(dets) < 7)
      dets <- cbind(dets, thickness) else
      dets[is.na(dets[,7]),7] <- thickness

    # should slumps be taken into account?
    if(length(slump) > 0)
      {
        d.adapt <- dets[,6]
        #d.lost <- c()
        d.lost <- NULL # has to be NULL since we don't know this var's final size
        for(i in 1:nrow(slump))
          {
            below.slump <- which(dets[,6] > max(slump[i,]))
            above.slump <- which(dets[,6] < min(slump[i,]))
            d.lost <- c(d.lost, which(!(1:nrow(dets) %in% c(above.slump, below.slump))))
            d.adapt[below.slump] <- d.adapt[below.slump] - (max(slump[i,])-min(slump[i,]))
          }
        dets[,6] <- d.adapt
        if(length(d.lost) > 0)
          dets <- dets[-d.lost,]
      }
    # check for common errors
    dets <- dets[,1:7]
    x <- 0
    for(i in 2:7) if(is.factor(dets[,i])) x <- 1
    if(x == 1)
      stop(paste("\n Some value fields in ", name, ".csv contain letters, please adapt", sep=""), call.=FALSE)
    if(length(dets[is.na(dets[,2]),2])+length(dets[is.na(dets[,3]),3]) != nrow(dets))
      stop(paste("\n Remove duplicate entries within the C14 and calendar fields in ", name, ".csv", sep=""), call.=FALSE)
    if(min(dets[,4]) <= 0)
      stop(paste("\n Errors of dates should be larger than zero. Please adapt ", name, ".csv", sep=""), call.=FALSE)
    dat$ID <- as.character(dets[,1])

    # correct for any reservoir effect
    dets[is.na(dets[,5]),5] <- 0
    dat$cage <- dets[,2] - dets[,5]
    dat$error <- dets[,4]

    # work in F14C for calibration
    dat$f.cage <- exp(-dat$cage/8033)
    dat$f.error <- dat$f.cage - exp(-(dat$cage+dat$error)/8033)

    # check if any 14C dates are (entirely or partly) beyond the calibration curve
    outside <- which(!is.na(dat$cage))
    rangecc <- c(min(calcurve[,2]-calcurve[,3]),max(calcurve[,2]+calcurve[,3]))
    outside <- outside[c(which(dat$cage[outside]-times*dat$error[outside] < rangecc[1]), which(dat$cage[outside]+times*dat$error[outside] > rangecc[2]))]
    if(length(outside) > 0)
      {
        truncate <- 0
        for(i in 1:length(outside)) # check if date lies only partly beyond the curve limits
          if((dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[1] &&
            dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[1]) ||
            (dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[2] &&
            dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[2]))
              truncate <- truncate + 1
        if(truncate > 0)
          cat("\n Warning, dates spanning beyond the calibration curve will be truncated! ")

        # remove dates which lie entirely outside the limits of the calibration curve
        outside <- outside[c(which(dat$cage[outside]+qnorm(1-(1-prob)/2)*dat$error[outside] < rangecc[1]), which(dat$cage[outside]-qnorm(1-(1-prob)/2)*dat$error[outside] > rangecc[2]))]
        if(length(outside) > 0)
          {
            cat("\n Warning, dates older than the calibration curve will be ignored! ")
            dets <- dets[-outside,]
            dat$cage <- dat$cage[-outside]
            dat$error <- dat$error[-outside]
            dat$f.cage <- dat$f.cage[-outside]
            dat$f.error <- dat$f.error[-outside]
            dat$outside <- dat$ID[outside]
            dat$ID <- dat$ID[-outside]
          }
      }

    # fill the 'dat' list with additional information
    dat$cal <- c(dets[,3], extradates)
    dat$res <- c(dets[,5], extradates)
    dat$depth <- c(dets[,6], extradates)
    dat$thick <- c(dets[,7], rep(thickness, length(extradates)))
    dat$BCAD <- BCAD

    # find distribution (calibrated if 14C) and point estimates for each date
    for(i in 1:length(dat$depth))
      {
        if(length(extradates) > 0 && i > nrow(dets))
          {
            tmp <- read.table(paste(dat$coredir, name, "_", extradates[i-nrow(dets)], ".txt", sep=""))
            calib <- cbind(tmp[,1], tmp[,2]/sum(tmp[,2]))
          } else
            if(is.na(dat$cage[[i]]))
              {
                age <- dat$cal[[i]]
                error <- dat$error[[i]]
                ageseq <- seq(age-(times*error), age+(times*error), by=yrsteps)
                calib <- cbind(ageseq, dnorm(ageseq, age, error))
              } else
                calib <- .caldist(dat$f.cage[[i]], dat$f.error[[i]], theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD)
        if(length(youngest) > 0) # truncate ages younger than a limit
          {
            if(BCAD) calib <- calib[which(calib[,1] <= youngest),] else
              calib <- calib[which(calib[,1] >= youngest),]
            if(length(calib) == 0)
              if(BCAD)
                calib <- cbind(seq(youngest-(3*yrsteps), youngest+yrsteps, length=5), c(0:3,0)/3) else
                calib <- cbind(seq(youngest-yrsteps, youngest+(3*yrsteps), length=5), c(0,3:0)/3)
          }
        dat$calib[[i]] <- calib  
        dat$hpd[[i]] <- .hpd(calib, prob=prob, hpdsteps=hpdsteps, yrsteps=yrsteps)
        dat$mid1[[i]] <- (dat$hpd[[i]][1] + dat$hpd[[i]][2*nrow(dat$hpd[[i]])])/2
        yrs <- calib[,1]
        dat$mid2[[i]] <- mean(c(max(yrs), min(yrs)))
        dat$wmn[[i]] <- weighted.mean(calib[,1], 1/calib[,2])
        dat$med[[i]] <- calib[max(which(cumsum(calib[,2]) <= .5)),1]
        dat$mode[[i]] <- calib[which(calib[,2] == max(calib[,2])),1][1]
      }

    if(storedat)
	  dets <<- dets
    dat
  }



# write files of the age-depth model, calibrated ranges, and settings
.write.clam <- function(dat, namedir,runname, calrange, name, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, cmyr, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, yrlab, dlab, plotrange, greyscale, chron, C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
  {
    # age-depth model; age estimates, accumulation rates and ranges for every analysed depth
    runnames <- c("_interpolated", "_polyn_regr", "_cubic_spline", "_smooth_spline", "_loess")
    calrange <- cbind(calrange, round(c(diff(calrange[,4])/diff(calrange[,1]), NA), decimals+2))
    if(cmyr)
      calrange[,5] <- 1/calrange[,5]
    calrange[,2:4] <- round(calrange[,2:4], decimals)
    ifelse(length(runname)==0, runname <- runnames[type], runname)
    if(depths.file && file.exists(dd <- paste(namedir, name, "/", name, "_depths.txt", sep="")))
      {
        dd <- read.table(dd)[,1]
        #this <- c()
		this <- numeric(length(dd))
        for(i in 1:length(dd))
          this[i] <- which(calrange[,1]==dd[i])[1] # find where the relevant ages are
	    write.table(calrange[this,], paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "acc.rate"), quote=FALSE, sep="\t")
      } else
       write.table(calrange, paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "accrate"), quote=FALSE, sep="\t")

    # calibrated ranges of all dates
    hpd.file <- file(paste(dat$coredir, name, "_calibrated.txt", sep=""), "w")
    cat(paste("Calibrated age ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
    for(i in 1:length(dat$depth))
      {
        cat(paste("\n\nDepth: ", dat$depth[[i]], "\nyrmin\tyrmax\tprobability\n"), file=hpd.file)
        hpds <- dat$hpd[[i]]
        for(j in 1:nrow(hpds))
          {
            for(k in 1:3) cat(hpds[j,k], "\t", file=hpd.file)
            cat("\n", file=hpd.file)
          }
      }
    close(hpd.file)

    # relevant settings and results
    set.file <- file(paste(dat$coredir, name, runnames[type], "_settings.txt", sep=""), "w")
    cat(paste("Settings (square brackets give names of the constants)\n\n",
      "Calibration curve: ", ccname,
      if(postbomb!=FALSE)
		paste(",", pbnames[postbomb], "for postbomb dates"),
      "\nAge-depth model: ",
      if(type==1) "linear interpolation between dated levels [type=1]" else
      if(type==2) ifelse(length(smooth)==0, "linear regression [type=2, smooth=c()]",
      paste("polynomial regression [type=2] of order", smooth, "[smooth]")) else
      if(type==3) "cubic spline [type=3]" else
      if(type==4) paste("smooth spline [type=4] with spar =", ifelse(length(smooth)<1, 0.3, smooth), "[smooth]") else
      if(type==5) paste("locally weighted spline [type=5] with span =", ifelse(length(smooth)<1, 0.75, smooth), "[smooth]"),
      if(wghts==1) "\nWeighted by the calibrated probabilities [wghts=1]",
      if(wghts==2) "\nWeighted by the errors (1/sdev^2) [wghts=2]",
      "\nCalculations at ", 100*prob, "% confidence ranges [prob=", prob, "]",
      "\nAmount of iterations: ", its, " [its]",
      "\nCalendar age point estimates for depths based on ",
      if(est==1) "weighted average of all age-depth curves [est=1]" else
      if(est==2) "midpoints of the hpd ranges of the age-depth curves [est=2]" else
      if(est==3) "midpoints of the hpd ranges of the dated levels [est=3]" else
      if(est==4) "weighted means of the dated levels [est=4]" else
      if(est==5) "medians of the dated levels [est=5]" else
      if(est==6) "modes/maxima/intercepts of the dated levels [est=6]",
      "\nCalendar scale used: ", if(BCAD) "cal BC/AD" else "cal BP",
      " [BCAD=", BCAD, "] at a resolution of ", yrsteps, " yr [yrsteps]",
      "\nAges were calculated every ", every, " [every] ", depth,
      " [depth], from ", min(depthseq), " [dmin] to ", max(depthseq), " [dmax] ", depth, sep=""), file=set.file)
      if(length(youngest) > 0) cat("\n\nDates with ages younger than", youngest, ifelse(BCAD, "BC/AD", "cal BP"), "were truncated", file=set.file)
      if(length(calibt)> 1) cat("\n\nInstead of assuming the standard Gaussian model, a student t distribution was used with t.a =", calibt[1], "and t.b =", calibt[2], "(see Christen and Perez 2009, Radiocarbon 51:1047-1059)", file=set.file)
    if(length(slump) == 2) cat("\n\nA slump was excised between", max(slump), "and", min(slump), depth, file=set.file)
    if(length(slump) > 2)
      {
        cat("\n\nSlumps were excised from ", file=set.file)
        sl <- array(sort(slump), dim=c(2, length(slump)/2))
        for(i in 1:ncol(sl))
          cat(sl[1,i], "to", sl[2,i], depth, if(i<ncol(sl)) "and ", file=set.file)
      }
    if(length(outliers) > 0)
      {
        cat("\n\nDates assumed outlying [outliers]: ", file=set.file)
        for(i in outliers) cat(i, " (", dat$ID[i], ") ", sep="", file=set.file)
      }
    if(length(ignore) > 0)
      {
        cat("\n\nDates ignored [ignore]: ", file=set.file)
        for(i in 1:length(ignore)) cat(ignore[i], " (", dat$ignore[i], ") ", sep="", file=set.file)
      }
    if(length(dat$outside) > 0)
      {
        cat("\n\nDates outside calibration curve and ignored: ", file=set.file)
        for(i in 1:length(dat$outside)) cat(dat$outside[i], " ", sep="", file=set.file)
      }
    cat(paste(
      if(length(hiatus) > 0)
        paste("\nA hiatus was inferred at", hiatus, depth, "[hiatus]"),
        "\n\nGoodness-of-fit (-log, lower is better): ", gfit,
      if(reversal) "\nSome age-depth reversals occurred"),
      if(remove.reverse) "\nAny models with age-depth reversals were removed",
      "\n\nProduced ", date(), sep="", file=set.file)
    close(set.file)

    if(plotpdf)
      {
        pdf(file=paste(dat$coredir, name, runname, ".pdf", sep=""))
        .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else NULL, C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
        dev.off()
      }
    if(plotpng)
      {
        png(filename = paste(dat$coredir, name, runname, ".png", sep=""))
        .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else NULL, C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
        dev.off()
      }
  }

  # If coredir is left empty, check for a folder named Cores in the current working directory, and if this doesn't exist, for a folder called clam_runs (make this folder if it doesn't exist yet and if the user agrees).
  # Check if we have write access. If not, tell the user to provide a different, writeable location for coredir. 
  assign_coredir <- function(coredir, core, ask=TRUE) {
    if(length(coredir) == 0) {
      if(dir.exists("Cores")) 
        coredir <- "Cores" else
          if(dir.exists("clam_runs"))
            coredir <- "clam_runs" else {
              coredir <- "clam_runs"
              ans <- readline(paste0("I will create a folder called ", coredir, ", is that OK? (y/n)  "))
              if(ask)
                if(tolower(substr(ans,1,1)) == "y")
                  wdir <- dir.create(coredir, FALSE) else
                    stop("No problem. Please provide an alternative folder location using coredir\n")
              if(!wdir)
                stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
          }    
    } else {
      if(!dir.exists(coredir))
          wdir <- dir.create(coredir, FALSE)
        if(!dir.exists(coredir)) # if it still doesn't exist, we probably don't have enough permissions
          stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
    }
    coredir <- .validateDirectoryName(coredir)
    cat("The run's files will be put in this folder: ", coredir, core, "\n", sep="")
    return(coredir)
  }


  # list the available cores
  clam_runs <- list.files("clam_runs/")

  .validateDirectoryName <- function(dir) {
    if(!dir.exists(dir))
      dir.create(dir, showWarnings=FALSE, recursive=TRUE)
    dir <- suppressWarnings(normalizePath(dir))
    lastchar <- substr(dir, nchar(dir), nchar(dir))
    if(lastchar != "/" & lastchar != "\\" & lastchar != "" & lastchar != "." )
      dir <- paste(dir, "/", sep="") # does this work in Windows?
    return(dir)
  }

  # function to load results into global environment
  # parameter position defaults to 1, which equals an assignment to the global environment
  .assign_to_global <- function(key, val, pos=1) {
    assign(key, val, envir=as.environment(pos) )
  }

