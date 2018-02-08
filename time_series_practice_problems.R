# Fish 507: Applied Time-Series Analysis in Fisheries and Environmental Sciences
# Week 1 & 2 Problems

# Load utility functions
source("/Volumes/workspaces/nate/public/coding_work/R/functions.R")
# Load packages
packages <- c("rstan", "RCurl", "xtable", "forecast", "stringr")
ipak(packages)

## You have been asked by a colleauge to help analyze some time series data she   ##
## collected as part of an experiment on the effects of light and nutrients on    ##
## the population dynamics of phytoplankton. Specifically, after controlling for  ##
## differences in light and temperature, she wants to know if the natural log of  ##
## population density can be modeled with some form of ARMA(p,q) model. The data  ## 
## are expressed as the number of cells per milliliter recorded every hour for    ##
## one week beginning at 8:00 AM on December 1, 2014; you can find them here:

## get phytoplankton data
pp <- "http://faculty.washington.edu/scheuerl/phytoDat.txt"
pDat <- read.table(pp)

## 1.1) Convert pDat, which is a data.frame object, into a ts object. This bit of  ##
## code might be useful to get you started:
## what day of 2014 is Dec 1st? Answer: it is the 335th day.
dBegin <- as.Date("2014-12-01")
dayOfYear <- (dBegin - as.Date("2014-01-01") + 1)
## Assign better column name
colnames(pDat) <- c("cells/ml")
## Convert to ts object. Arguments of ts handle the date/time wrangling.
pDat <- ts(data=pDat, frequency=24,
          start=c(dayOfYear, 8))

## 1.2) Plot the time series of phytoplankton density and provide a brief description ##
## of any notable features:

## plot the ts
plot.ts(pDat, xlab="Time (day of year)", ylab="Cells/ml")
## NOTE: Data appear to be non-stationary in that they have obvious seasonal signal   ##
## corresponding to each 24-hr period. Data also appear to be non-Gaussian (i.e. they ##
## occur on the positive, real interval).                                             ##

## 1.3) Although you do not have the actual measurements for the specific temperature  ##
## and light regimes used in the experiment, you have been informed that they follow  ##
## a regular light/dark period with accompanying warm/cool temperatures. Thus,        ##
## estimating a fixed seasonal effect is justifiable. Also, the instrumentation is    ##
## precise enough to preclude any systematic change in measurements over time (i.e.,  ##
## you can assume mt = 0 for all t). Obtain the time series of the estimated          ##
## log-density of phytoplankton absent any hourly effects caused by variation in      ##
## temperature or light. (Hint: You will need to do some decomposition.)              ##

## Log transform phytoplankton density:
log.pDat <- log(pDat)
## Need to average the estimates of seasonal effect for each period over the whole study:
## length of ts:
len <- length(log.pDat)
## frequency of ts:
freq <- frequency(log.pDat)
## number of periods
periods <- len %/% freq
## index of cumulative period
index <- seq(1, len, by=freq) - 1
## get mean by period:
mm <- numeric(freq)
for (i in 1:freq) {
  mm[i] <- mean(log.pDat[index+i], na.rm = TRUE)
}
## subtract mean/period to make overall mean=0
mm <- mm - mean(mm)
## plot the daily seasonal effects
plot.ts(mm, ylab="Seasonal effect", xlab="Day", cex=1)
## create entire time series of seasonal effects:
## create ts object for season:
log.pDat.seas <- ts(rep(mm, periods+1)[seq(len)], 
                        start=start(log.pDat), 
                        frequency=freq)
## random errors over time:
log.pDat.err <- log.pDat - log.pDat.seas
## plot the observation ts & seasonal effect ts
plot(cbind(log.pDat,log.pDat.seas,log.pDat.err),main="",yax.flip=TRUE)

## Alternative method: use decompose
log.pDat.decomp <- decompose(log.pDat)
## plot the obs ts, trend & seasonal effect
plot(log.pDat.decomp, yax.flip=TRUE)

## 1.4) Use diagnostic tools to identify the possible order(s) of ARMA model(s) that  ##
## most likely describes the log of population density for this particular experiment.##
## Note that at this point you should be focusing your analysis on the results        ##
## obtained in Question 3.

## calculate ACF:
## NOTE: lag.max is set to the total number of samples
log.pDat.acf <- acf(log.pDat, lag.max = 168)
## Alternative custom ACF plot function:
plot.acf <- function(ACFobj) {
  rr <- ACFobj$acf[-1]
  kk <- length(rr)
  nn <- ACFobj$n.used
  plot(seq(kk),rr,type="h",lwd=2,yaxs="i",xaxs="i",
       ylim=c(floor(min(rr)),1),xlim=c(0,kk+1),
       xlab="Lag",ylab="Correlation",las=1)
  abline(h=-1/nn+c(-2,2)/sqrt(nn),lty="dashed",col="blue")
  abline(h=0)
}
## Pass ACF object to custom plot function:
plot.acf(log.pDat.acf)

## calculate PACF:
log.pDat.pacf <- pacf(log.pDat, lag.max = 168)
## Alternative custom ACF plot function:
plot.pacf <- function(PACFobj) {
  rr <- PACFobj$acf
  kk <- length(rr)
  nn <- PACFobj$n.used
  plot(seq(kk),rr,type="h",lwd=2,yaxs="i",xaxs="i",
       ylim=c(floor(min(rr)),1),xlim=c(0,kk+1),
       xlab="Lag",ylab="PACF",las=1)
  abline(h=-1/nn+c(-2,2)/sqrt(nn),lty="dashed",col="blue")
  abline(h=0)
}
## Pass PACF object to custom plot function:
plot.pacf(log.pDat.pacf)
## Use arima.sim function to simulate data:
## NOTE: model argument of arima.sim is list of following elements:
## order a vector of length 3 containing the ARIMA(p,d,q) order ar a vector of length p containing the AR(p) coefficients ##
## ma a vector of length q containing the MA(q) coefficients                                                              ##
## sd a scalar indicating the std dev of the Gaussian errors                                                              ##
## NOTE: you can pass arima.sim your own time series of random errors                                                     ##

## Use arima function to fit an ARIMA model to a univariate time series
## IMPORTANT ARGS:
## x a univariate time series
## order a vector of length 3 specifying the order of ARIMA(p,d,q) model

## Use auto.arima function to search over all possible orders of ARIMA models you specify  to see which fits most closely
## find best ARMA(p,q) model
auto.arima(log.pDat.err, max.p=3, max.q=3, seasonal=FALSE, trace = 1)
## “best” model is an AR(1) model that includes a non-zero intercept with following equation:
## x_t = 1.15 + 0.56*(x_{t-1} - 1.15) + w_t with w_t ~ N(0,0.085)