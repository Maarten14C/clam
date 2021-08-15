

#' @name deptime.depth
#' @title  Calculates *for each iteration* the slope of a straight curve between depths
#'  just above and below the desired point.
#' @description Calculates *for each iteration* the slope of a straight curve between depths
#'  above and below the desired point. Requires sufficiently dense density of depths, e.g. \code{yrsteps=1}.
#' @details 
#' To calculate sedimentation times at a depth. Before running this, run your core in clam and store the data, 
#' so, make sure to set \code{storedat=TRUE}. 
#' Renamed from previous accrate.depth function to avoid confusion with accrate.depth function of rbacon.
#' @param depth The depth for which accumulation rate estimates should be calculated.
#' @param yrcm Calculate in years per cm, or alternatively in cm per yr.
#' @param prob  Probability level at which to calculate the ranges.
#' @author Maarten Blaauw
#' @return Returns (invisibly) the modelled deposition times for a specific depths, a histogram and confidence ranges.
#' @examples
#'   clam(coredir=tempdir(), storedat=TRUE) 
#'   dp <- deptime.depth(20)
#'   summary(dp)
#'   deptime.depth(20, FALSE) # to calculate accumulation rates in cm/yr
#' 
#' @export
deptime.depth <- function(depth, yrcm=TRUE, prob=.95) {
  chron <- get('chron')
  calrange <- get('calrange')
  if(depth <= min(calrange[,1]) || depth >= max(calrange))
    stop("Deposition times cannot be calculated for the top or bottom of the core. Please check the manual", call.=FALSE)
  d <- max(which(calrange[,1] <= depth))
  if(yrcm)
    accrate <- (chron[d+1,]-chron[d-1,]) / (calrange[d+1,1]-calrange[d-1,1]) else
      accrate <- (calrange[d+1,1]-calrange[d-1,1]) / (chron[d+1,]-chron[d-1,])
  acc <- density(accrate)
  plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
  abline(h=0)
  o <- order(acc$y, decreasing=TRUE)
  acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
  acc <- range(acc[acc[,2] <= prob,1])
  rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
  cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")
  invisible(accrate)
}



#' @name deptime.age 
#' @title Calculates the slope of a straight curve at the desired age.
#' @description Calculates *for each iteration* the slope of a straight curve between 
#' depths above and below the desired age. Requires sufficiently dense density of depths, e.g. \code{steps=1}.
#' @details 
#' To calculate deposition times at an age. Before doing this, run your core in clam and store the data, 
#' so, make sure the option \code{storedat=TRUE}.
#' Renamed from previous accrate.age function to avoid confusion with accrate.age function of rbacon.
#' @param age Age to calculate deposition time (years per cm).
#' @param yrcm Calculate in years per cm, or alternatively in cm per yr.
#' @param prob Probability level at which to calculate the ranges.
#' @author Maarten Blaauw
#' @return Returns (invisibly) the modelled deposition times at a specific age, a histogram and confidence ranges.
#' @examples 
#'   clam(coredir=tempdir(), storedat=TRUE)
#'   dp <- deptime.age(5000)
#'   summary(dp)
#'   deptime.age(5000, yrcm=FALSE) # to calculate sedimentation times in cm/yr, so accumulation rates
#' @export
deptime.age <- function(age, yrcm=TRUE, prob=.95) {
  chron <- get('chron')
  calrange <- get('calrange')

  accrate <- numeric(ncol(chron))
  # accrate <- c()
  for(i in 1:ncol(chron)) {
    a <- max(which(chron[,i] <= age))
    if(yrcm)
      accrate <- c(accrate, (chron[a+1,i]-chron[a-1,i]) / (calrange[a+1,1]-calrange[a-1,1])) else
        accrate <- c( accrate, (calrange[a+1,1]-calrange[a-1,1]) / (chron[a+1,i]-chron[a-1,i]))
  }
  acc <- density(accrate)
  plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
  abline(h=0)
  o <- order(acc$y, decreasing=TRUE)
  acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
  acc <- range(acc[acc[,2] <= prob,1])
  rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
  cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")
  invisible(as.double(accrate))
}



#' @name plot_proxies 
#' @title Produce a plot of proxy values against calendar age.
#' @description  Produce a plot of proxy values against calendar age.
#' @details 
#' Only works after running clam on the core using \code{proxies=TRUE}. Requires a file containing the core depths as the first column, 
#' and any proxy values on subsequent columns. Values should be separated by comma's. The file should be stored as a .csv file in the core's directory.
#' @param prox  Position of the proxy that should be plotted, e.g. \code{1} for the first proxy in the file.
#' @param errors Plot an error envelope.
#' @param proxcol Colour of the error envelope.
#' @param revyr Direction of the calendar scale (\code{revyr=TRUE} will reverse the calendar scale from the default \code{FALSE}).
#' @author Maarten Blaauw
#' @return A plot of the age model function with proxies.
#' @examples 
#' clam(coredir=tempdir(), proxies=TRUE)
#' plot_proxies(3)
#' plot_proxies(3, revyr=FALSE)
#' @export
plot_proxies <- function(prox, errors=TRUE, proxcol=grey(0.5), revyr=TRUE) {
  dat <- get('dat') #JEV warning
  calrange <- get('calrange') #JEV warning
  
  prx <- dat$proxies
  if(length(prox)>1) layout(matrix(1:length(prox), ncol=1))
  for(j in 1:length(prox)) {
    pr <- prx[which(!is.na(prx[,prox+1])),]
    ages <- array(0, dim=c(nrow(pr),3))
    for(i in 1:nrow(pr))
      ages[i,] <- calrange[which(calrange[,1]==pr[i,1]),c(2,3,4)]
    xlim <- range(ages)
    if(!dat$BCAD) xlim <- rev(xlim)
    if(revyr) xlim <- rev(xlim)
    plot(ages[,3], pr[,prox+1], type="n", xlim=xlim, xlab=ifelse(dat$BCAD, "cal BC/AD", "cal BP"), ylab=names(pr)[prox+1])
    if(errors)
      for(i in 2:nrow(pr))
        polygon(c(ages[(i-1):i,1], ages[i:(i-1),2]), c(pr[c((i-1):i, i:(i-1)),prox+1]), col=proxcol, border=proxcol)
    lines(ages[,3], pr[,prox+1])
  }
  layout(1)
}
