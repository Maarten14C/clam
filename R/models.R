# sample point age estimates from the calibrated distributions ('its' times)
# the probability of a year being sampled is proportional to its calibrated probability
.smpl <- function(its, depths, calibs, Est) {
  smp <- array(1, dim=c(length(depths), 1+its, 2))
  #cat("\nEst=", length(Est), dim(smp), "\n")
  smp[,1,1] <- Est
  for(i in 1:length(calibs)) {
    thiscalib <- as.matrix(calibs[[i]], ncol=2) 
    sampled <- sample(1:nrow(thiscalib), its, prob=thiscalib[,2], TRUE)
    smp[i,(1:its)+1,] <- thiscalib[sampled,]
  }
  return(smp)
}

# akin to Heegaard et al.'s mixed effect modelling, but using calibrated dates
.mixed.effect <- function(its, depths, cals, cages, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) {
    cat("\n Mixed effect modelling, this will take some time")
    smp <- array(1, dim=c(length(depths), 1+its, 2))
    smp[,1,1] <- Est
    for(i in 1:length(cals))
      if(!is.na(cals[i])) {
        if(length(calibt) == 0) {
          x <- rnorm(its, cals[i], errors[i])
          smp[i,(1:its)+1,] <- c(x, dnorm(x, cals[i], errors[i]))
        } else {
            x <- (cals[i]-10*errors[i]) : (cals[i]+10*errors[i])
      #      x <- cbind(x, .calibt(calibt[1], calibt[2], cals[i], errors[i], x, 0))
            x <- cbind(x, dnorm(x, cals[i], errors[i])) # since calibt is not working
            o <- order(x[,2], decreasing=TRUE)
            x <- cbind(x[o,1], cumsum(x[o,2])/sum(x[,2]))
            sampled.x <- max(which(x[,2] <= runif(1, 0, max(x[,2]))))
            smp[i,(1:its)+1,] <- x[sampled.x,]
          }
      } else
         for(j in 1:its) {
           if(j/(its/3) == round(j/(its/3))) cat(".")
           yr <- rnorm(1, cages[i], errors[i])
           f.yr <- exp(-yr/8033)
           f.error <- f.yr - exp(-(yr+errors[i])/8033)
           yr <- cbind(theta, dnorm(f.mu, f.yr, sqrt(f.error^2+f.sigma^2)))
           yr <- yr[yr[,2]>0,]
           yr <- approx(yr[,1], yr[,2], seq(min(yr[,1]), max(yr[,1]), by=yrsteps))
           smp.yr <- sample(length(yr$x), 1, prob=yr$y)
           smp[i,j+1,] <- c(yr$x[smp.yr], yr$y[smp.yr])
         }
  return(smp)
}


# interpolate linearly between the data (default)
.interp <- function(depthseq, depths, its, chron, smp) {
  cat(" Interpolating, sampling")
  for(i in 1:its) {
    temp <- approx(depths, smp[,i,1], depthseq, ties=mean)$y

    # allow for extrapolation... dangerous!
    if(min(depthseq) < min(depths)) {
      minus <- which(depthseq < min(depths))
      slope <- diff(temp)[max(minus)+1]/diff(depthseq)[max(minus)+1]
      temp[minus] <- temp[max(minus)+1] + slope * (depthseq[minus] - min(depths))
    }
    if(max(depthseq) > max(depths)) {
      maxim <- which(depthseq > max(depths))
      slope <- diff(temp)[min(maxim)-2]/diff(depthseq)[min(maxim)-2]
      temp[maxim] <- temp[min(maxim)-1] + slope * (depthseq[maxim] - max(depths))
    }
    chron[,i] <- temp
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  return(chron)
}


# polynomial regressions of certain order through the data (default linear, y=ax+b)
.poly <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp) {
  if(length(smooth)==0)
    message(" Using linear regression, sampling") else
      message(paste0(" Using polynomial regression (degree ", smooth, "), sampling"))
  if(wghts==0) w <- NULL else w <- 1/errors^2
  for(i in 1:its) {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(lm(smp[,i,1] ~ poly(depths, max(1, smooth)), weights=w), data.frame(depths=depthseq))
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  return(chron)
}


# fit cubic spline interpolations through the data
.spline <- function(depthseq, smooth, depths, its, chron, smp) {
  if(length(smooth) < 1) smooth <- .3
  message(" Using cubic spline sampling")
  for(i in 1:its) {
    chron[,i] <- spline(depths, smp[,i,1], xout=depthseq)$y
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  return(chron)
}


# fit cubic smoothed splines through the data, with smoothing factor
.smooth <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp) {
  if(length(smooth) < 1) smooth <- .3
  message(paste0(" Using smoothing spline (smoothing ", smooth, "), sampling"))
  #if(wghts==0) w <- c() else w <- 1/errors^2
  if(wghts==0) w <- NULL else w <- 1/errors^2
  for(i in 1:its) {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(smooth.spline(depths, smp[,i,1], w=w, spar=smooth), depthseq)$y
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  return(chron)
}


# fit locally weighted (1/errors^2) splines through the data, with smoothing factor
.loess <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp) {
  if(length(smooth) < 1) smooth <- .75
  message(paste0(" Using loess (smoothing ", smooth, "), sampling"))
  #if(wghts==0) w <- c() else w <- 1/errors^2
  if(wghts==0) w <- NULL else w <- 1/errors^2
  for(i in 1:its) {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(loess(smp[,i,1] ~ depths, weights=w, span=smooth), depthseq)
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  return(chron)
}

# calculate goodness-of-fit (small number, so calculate its -log)
.gfit <- function(theta, f.mu, f.sigma, dat, calrange, outliers) {
  # gfit <- c()
  if(length(outliers) > 0) {
    dat$cage <- dat$cage[-outliers]
    dat$error <- dat$error[-outliers]
    dat$cal <- dat$cal[-outliers]
    dat$model <- dat$model[-outliers]
  }
  gfit <- pnorm(dat$cal, dat$model, dat$error^2)
  if(length(c14 <- which(!is.na(dat$cage))) > 0) { # if there are radiocarbon dates
    gfit.c <- approx(theta, f.mu, dat$model[c14])$y # C14 age at cc of modelled cal date
    f.cage <- exp(-dat$cage[c14]/8033)
    f.error <- exp(-(dat$cage[c14]-dat$error[c14])/8033) - f.cage
    gfit.var <- f.error^2 + approx(theta, f.sigma, dat$model[c14])$y^2
    gfit[c14] <- pnorm(f.cage, gfit.c, sqrt(gfit.var)) # deviation between measured and cc ages
  }
  dat$gfit <- -sum(log(gfit[!is.na(gfit)]))
}

