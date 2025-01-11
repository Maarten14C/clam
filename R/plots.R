
.ageplot <- function(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, yrlab, dlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, chron, C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty="l", mar, mgp, ash=FALSE) {
    # set up initial parameters
    if(length(dlab)==0) dlab <- paste0("Depth (", depth, ")")
    ifelse(BCAD || !revyr, yr.lim <- c(yrmin, yrmax), yr.lim <- c(yrmax, yrmin))
    if(revd) d.lim <- c(dmax, dmin) else d.lim <- c(dmin, dmax)

    par(xaxt="s", xaxs="r", yaxt="s", yaxs="r", bty=bty, mar=mar, mgp=mgp, font=2)
    if(revaxes) plot(0, type="n", ylim=yr.lim, xlim=d.lim, xlab=dlab, ylab=yrlab) else
      plot(0, type="n", xlim=yr.lim, ylim=d.lim, xlab=yrlab, ylab=dlab)
    if(plotname) legend("topleft", name, bty="n")

    # draw histograms of all age-depth models. Off by default, time-consuming!
    if(length(greyscale)==1) {
      plotrange <- FALSE
      depgr=seq(dmin, dmax, length=greyscale)
      for(i in 2:greyscale) {
        temp <- density(chron[max(which(calrange[,1]<=depgr[i])),2:ncol(chron)], n=greyscale)
        if(revaxes)
          image(c(depgr[i-1], depgr[i]), temp$x, matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE) else
            image(temp$x, c(depgr[i-1], depgr[i]), matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE)
      }
    }

    # draw the age-depth models, per section if hiatuses were inferred
    if(length(hiatus) > 0) {
      if(length(slump) == 0)
        hiatusseq <- sort(c(range(depthseq), hiatus)) else
          hiatusseq <- sort(c(range(depthseq, depthseq+sum(slump[,2]-slump[,1])), hiatus))
      for(i in 2:length(hiatusseq)) {
        sec <- calrange[min(which(calrange[,1] > hiatusseq[i-1])):max(which(calrange[,1] < hiatusseq[i])),]
        pol <- cbind(c(sec[,2], rev(sec[,3])), c(sec[,1], rev(sec[,1])))
        if(plotrange)
          if(revaxes)
            polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
              polygon(pol, col=rangecol, border=rangecol)
        if(revaxes)
          lines(sec[,1], sec[,4], lwd=2, col=bestcol) else
            lines(sec[,4], sec[,1], lwd=2, col=bestcol)
        if(revaxes)
          abline(v=hiatus, col="grey", lty="dashed") else
            abline(h=hiatus, col="grey", lty="dashed")
      }
    } else {
        pol <- cbind(c(calrange[,2], rev(calrange[,3])), c(calrange[,1], rev(calrange[,1])))
        if(plotrange)
          if(revaxes)
            polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
              polygon(pol, col=rangecol, border=rangecol)
        if(revaxes)
          lines(calrange[,1], calrange[,4], lwd=2, col=bestcol) else
            lines(calrange[,4], calrange[,1], lwd=2, col=bestcol)
      }

    # draw slumps if these were given
    if(length(slump) > 0)
      for(i in 1:nrow(slump))
        if(revaxes)
          rect(min(slump[i,]), min(yr.lim)-1e4, max(slump[i,]), max(yr.lim)+1e4, col=slumpcol, border=slumpcol) else
            rect(min(yr.lim)-1e4, min(slump[i,]), max(yr.lim)+1e4, max(slump[i,]), col=slumpcol, border=slumpcol)

    # draw the calibrated distributions of the dates
    top <- 1
    for(i in 1:length(dat$depth))
      top <- min(top, max(dat$calib[[i]][,2])) # find the lowest peak

    if(calhght > 0)
      for(i in 1:length(dat$depth)) {
        if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
        pol <- dat$calib[[i]] # already normalised to 1
        if(ash) pol[,2] <- pol[,2]/max(pol[,2])/1e3 # draw all same height
        pol[pol[,2] > maxhght,2] <- maxhght
        pol[,2] <- calhght*(dmax-dmin)*pol[,2]/(top*100)
        pol <- cbind(c(pol[,1], rev(pol[,1])),
          c(dat$depth[[i]]-pol[,2], dat$depth[[i]]+mirror*rev(pol[,2])))
        if(revaxes) polygon(pol[,2], pol[,1], col=col, border=col) else
          polygon(pol, col=col, border=col)
      }

    # draw the calibrated ranges of the dates
    for(i in 1:length(dat$depth)) {
      if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
      for(j in 1:nrow(dat$hpd[[i]]))
        if(revaxes)
          rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], lwd=1, lend=2, col=col, border=NA) else
            rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, lwd=1, lend=2, col=col, border=NA)
    }
    if(length(outliers) >0 ) {
      for(i in outliers)
        for(j in 1:nrow(dat$hpd[[i]]))
          if(revaxes)
            rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], col=outcol, border=outcol, lwd=1, lend=2) else
              rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, col=outcol, border=outcol, lwd=1, lend=2)
        if(revaxes)
          points(dat$depth[outliers], dat$mid1[outliers], cex=outlsize, pch=4, col=outcol) else
            points(dat$mid1[outliers], dat$depth[outliers], cex=outlsize, pch=4, col=outcol)
    }
  }



#' @name add.dates
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets.
#' @param mn Reported mean of the date. Can be multiple dates.
#' @param sdev Reported error of the date. Can be multiple dates.
#' @param depth Depth of the date.
#' @param cc The calibration curve to use: \code{cc=1} for IntCal20 (northern hemisphere terrestrial), \code{cc=2} for Marine20 (marine), \code{cc=3} for SHcal20 (southern hemisphere terrestrial), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param above Threshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param exx Exaggeration of probability distribution plots. Defaults to \code{exx=50}.
#' @param normal By default, Bacon uses the student's t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param normalise By default, the date is normalised to an area of 1 (\code{normalise=TRUE}).
#' @param t.a The dates are treated using the student's t distribution by default (\code{normal=FALSE}).
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param t.b The dates are treated using the student's t distribution by default (\code{normal=FALSE}).
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010).
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file).
#' For symmetry reasons, t.a must always be equal to t.b-1.
#' @param age.res Resolution of the date's distribution. Defaults to \code{date.res=100}.
#' @param times The extent of the range to be calculated for each date. Defaults to \code{times=20}.
#' @param col The colour of the ranges of the date. Default is semi-transparent red: \code{col=rgb(1,0,0,.5)}.
#' @param border The colours of the borders of the date. Default is semi-transparent red: \code{border=rgb(1,0,0,0.5)}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}.
#' @author Maarten Blaauw, J. Andres Christen
#' @return A date's distribution, added to an age-depth plot.
#' @examples
#'   base_temp_dir <- tempdir()
#'   clam_dir <- file.path(base_temp_dir, "clam_runs")
#'   dir.create(clam_dir, recursive = TRUE, showWarnings = FALSE)
#'   clam(, coredir=clam_dir, ask=FALSE)
#'   add.dates(5000, 100, 60)
#' @export
add.dates <- function(mn, sdev, depth, cc=1, above=1e-3, exx=50, normal=TRUE, normalise=TRUE, t.a=3, t.b=4, age.res=100, times=20, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE)  {
  if(cc > 0)
    cc <- rintcal::ccurve(cc) # was ccurve(cc) Sep 2024
  
  for(i in 1:length(mn)) {
    yrs <- seq(mn[i]-times*sdev[i], mn[i]+times*sdev[i], length=age.res)
    if(length(cc) < 2)
      cc <- cbind(yrs, yrs, rep(0, length(yrs)))
  ages <- approx(cc[,1], cc[,2], yrs)$y
  errors <- approx(cc[,1], cc[,3], yrs)$y

  if(normal)
    probs <- dnorm(ages, mn[i], sqrt(sdev[i] + errors)^2) else
      probs <- (t.b + (mn[i]-ages)^2  / (2*(sdev[i]^2 + errors^2))) ^ (-1*(t.a+0.5))
  if(normalise)
    probs <- probs / sum(probs)
  these <- which(probs >= above)
  if(length(these) > 0) {
    yrs <- yrs[these]
    probs <- probs[these]
  }

  if(!up)
    up <- -1
  if(BCAD)
    yrs <- 1950 - yrs
  if(mirror)
    pol <- cbind(c(yrs, rev(yrs)), depth[i] + exx*c(probs, -rev(probs))) else
      pol <- cbind(c(min(yrs), yrs, max(yrs)), depth[i] - up*exx*c(0, probs,  0))
  if(rotate.axes)
    pol <- pol[,2:1]
  polygon(pol, col=col, border=border)
  }
}


