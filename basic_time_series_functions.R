# Basic Time Series Functions In R
# Fish 507: Applied Time-Series Analysis in Fisheries and Environmental Sciences

# Load utility functions
source("/Volumes/workspaces/nate/public/coding_work/R/functions.R")
# Load packages
packages <- c("rstan", "RCurl", "xtable", "forecast", "stringr")
ipak(packages)

## TIME SERIES PLOTS ##
## get CO2 data from Mauna Loa observatory
ww1 <- "ftp://aftp.cmdl.noaa.gov/products/"
ww2 <- "trends/co2/co2_mm_mlo.txt"
CO2 <- read.table(text=getURL(paste0(ww1,ww2)))[,c(1,2,5)]
## assign better column names
colnames(CO2) <- c("year","month","ppm")

## create a time series (ts) object from the CO2 data
co2 <- ts(data=CO2$ppm, frequency=12,
          start=c(CO2[1,"year"],CO2[1,"month"]))
## plot the ts
plot.ts(co2, ylab=expression(paste("CO"[2]," (ppm)")))

## get N Hemisphere land & ocean temperature anomalies from NOAA
ww1 <- "https://www.ncdc.noaa.gov/cag/time-series/"
ww2 <- "global/nhem/land_ocean/p12/12/1880-2014.csv"
Temp <- read.csv(text=getURL(paste0(ww1,ww2)), skip=4)
## create ts object
tmp <- ts(data=Temp$Value, frequency=12, start=c(1880,1))

## intersection (only overlapping times)
datI <- ts.intersect(co2,tmp)
## dimensions of common-time data
dim(datI)

## union (all times)
datU <- ts.union(co2,tmp)
## dimensions of all-time data
dim(datU)

## plot the ts
plot(datI, main="", yax.flip=TRUE)

## DECOMPOSITION OF TIME SERIES             ##
## 1) trend (t)                             ##
## 2) seasonal effects (mt)                 ##
## 3) random errors (et)                    ##
## simple additive decomposition model: xt  ##
## xt = mt +st +et                          ##

## weights for moving avg
fltr <- c(1/2,rep(1,times=11),1/2)/12

## estimate of trend
co2.trend <- filter(co2, filter=fltr, method="convo", sides=2)
## plot the trend
plot.ts(co2.trend, ylab="Trend", cex=1)

## seasonal effect over time
co2.1T <- co2 - co2.trend

## plot the monthly seasonal effects
plot.ts(co2.1T, ylab="Seasonal effect", xlab="Month", cex=1)

## length of ts
ll <- length(co2.1T)
## frequency (ie, 12)
ff <- frequency(co2.1T)
## number of periods (years); %/% is integer division
periods <- ll %/% ff
## index of cumulative month
index <- seq(1,ll,by=ff) - 1
## get mean by month
mm <- numeric(ff)

for(i in 1:ff) {
  mm[i] <- mean(co2.1T[index+i], na.rm=TRUE)
}
## subtract mean to make overall mean=0
mm <- mm - mean(mm)

## plot the monthly seasonal effects
plot.ts(mm, ylab="Seasonal effect", xlab="Month", cex=1)

## create ts object for season
co2.seas <- ts(rep(mm, periods+1)[seq(ll)],
               start=start(co2.1T),
               frequency=ff)

## random errors over time
co2.err <- co2 - co2.trend - co2.seas

## plot the obs ts, trend & seasonal effect
plot(cbind(co2,co2.trend,co2.seas,co2.err),main="",yax.flip=TRUE)

## decomposition of CO2 data
co2.decomp <- decompose(co2)

## plot the obs ts, trend & seasonal effect
plot(co2.decomp, yax.flip=TRUE)

## DIFFERENCING TO REMOVE A TREND OR SEASONAL EFFECTS ##
## (alternative to decomposition)                     ##
## Difference operator: ∇xt = xt −xt−1                ##

## twice-difference the CO2 data
co2.D2 <- diff(co2, differences=2)
## plot the differenced data
plot(co2.D2, ylab=expression(paste(nabla^2,"CO"[2])))

## Note: trend has been successfully removed, but seasonal effect is still obvious. Therefore, difference
## this already-differenced series at lag-12 because data was collected monthly:
co2.D2D12 <- diff(co2.D2, lag=12)
## plot the newly differenced data
plot(co2.D2D12,
     ylab=expression(paste(nabla,"(",nabla^2,"CO"[2],")")))
## Note: this operation successfully produces a time series that looks like random errors without any 
## obvious trend or seasonal components.

## CORRELATION WITHIN AND AMONG TIME SERIES ##
## Examine the correlation structure of the original data or random errors from a decomposition model ##
## to help us identify possible form(s) of (non)stationary model(s) for the stochastic process.       ##
## AUTOCORRELATION FUNCTION (ACF) ##

## correlogram of the CO2 data
acf(co2, lag.max=36)

## As an alternative to the plotting utility in acf, let’s define a new plot function for acf objects ##
## with some better features:                                                                         ##
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

## Assign the result of acf to a variable and then use the infor- mation contained therein to plot  ##
## the correlogram with our new plot function.                                                      ##
## acf of the CO2 data
co2.acf <- acf(co2, lag.max=36)
## correlogram of the CO2 data
plot.acf(co2.acf)

## Look at the ACF for some deter- ministic time series, which will help you identify interesting   ##
## properties (e.g.,trends, seasonal effects) in a stochastic time series, and account for them in  ##
## time series models.                                                                              ##
## Example: straight line:
## length of ts
nn <- 100
## create straight line
tt <- seq(nn)
## set up plot area
par(mfrow=c(1,2))
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## get ACF
line.acf <- acf(tt, plot=FALSE)
## plot ACF
plot.acf(line.acf)

## Example: sine wave:
## create sine wave
tt <- sin(2*pi*seq(nn)/12)
## set up plot area
par(mfrow=c(1,2))
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## get ACF
sine.acf <- acf(tt, plot=FALSE)
## plot ACF
plot.acf(sine.acf)

## Example: sine wave with a linear downward trend:
## create sine wave with trend
tt <- sin(2*pi*seq(nn)/12) - seq(nn)/50
## set up plot area
par(mfrow=c(1,2))
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## get ACF
sili.acf <- acf(tt, plot=FALSE)
## plot ACF
plot.acf(sili.acf)

## PARTIAL AUTOCORRELATION FUNCTION (PACF) ##
## PACF of the CO2 data
pacf(co2, lag.max=36)

## Useful custom PACF plotting function: 
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

## Try the custom PACF plot function:
## PACF of the CO2 data
co2.pacf <- pacf(co2)
## correlogram of the CO2 data
plot.acf(co2.pacf)

## CROSS-CORRELATION FUNCTION (CCF) ##
## NOTE: pay particular attention to which variable you call y (i.e.,the “response”) and which 
## you call x (i.e., the “predictor”).

## Example: sunspot activity vs number of lynx trapped in Canada (Moran, 1949)
## get the matching years of sunspot data
suns <- ts.intersect(lynx,sunspot.year)[,"sunspot.year"]
## get the matching lynx data
lynx <- ts.intersect(lynx,sunspot.year)[,"lynx"]
## plot time series
plot(cbind(suns,lynx), yax.flip=TRUE)
## CCF of sunspots and lynx
ccf(suns, log(lynx), ylab="Cross-correlation")

## WHITE NOISE (WN) ##
## A time series {wt} is a discrete white noise series (DWN) if the w1,w1,...,wt are independent ##
## and identically distributed (IID) with a mean of zero. For most of the examples in this       ##
## course we will assume that the wt ∼ N(0, q), and therefore we refer to the time series {wt}   ##
## as Gaussian white noise.                                                                      ##
## SIMULATE WN: 
set.seed(123)
## random normal variates
GWN <- rnorm(n=1000, mean=5, sd=0.2)
## random Poisson variates
PWN <- rpois(n=1000, lambda=20)
## set up plot region
par(mfrow=c(1,2))
## plot normal variates with mean
plot.ts(GWN)
abline(h=5, col="blue", lty="dashed")
## plot Poisson variates with mean
plot.ts(PWN)
abline(h=20, col="blue", lty="dashed")

## Examine the ACF for the 2 white noise series and see if there is, in fact, zero autocorrelation ##
## for lags ≥ 1.                                                                                   ##
## set up plot region
par(mfrow=c(1,2))
## plot normal variates with mean
acf(GWN, main="", lag.max=20)
## plot Poisson variates with mean
acf(PWN, main="", lag.max=20)

## RANDOM WALKS (RW) ##
## Random walks are the most simple non-stationary time series model. A random walk is a time     ##
## series {xt} where xt =xt−1+wt, and wt is a discrete white noise series where all values are    ##
## independent and identically distributed (IID) with a mean of zero. In practice, we will almost ##
## always assume that the wt are Gaussian white noise, such that wt ∼ N(0, q). We will see later  ##
## that a random walk is a special case of an autoregressive model.                               ##
## SIMULATE A RW:
## set random number seed
set.seed(123)
## length of time series
TT <- 100
## initialize {x_t} and {w_t}
xx <- ww <- rnorm(n=TT, mean=0, sd=1)
## compute values 2 thru TT
for(t in 2:TT) { xx[t] <- xx[t-1] + ww[t] }
## Plot simulated ts and its ACF:
## setup plot area
par(mfrow=c(1,2))
## plot line
plot.ts(xx, ylab=expression(italic(x[t])))
## plot ACF
plot.acf(acf(xx, plot=FALSE))

## ALTERNATIVE METHOD FOR SIMULATING RW:
## Based on fact that value of RW process at time step t is sum of all random errors up through time t. 
## Use cumulative summation of vector x over its entire length:
## simulate RW
x2 <- cumsum(ww)
## plot ts to check if alternative method worked.
## setup plot area
par(mfrow=c(1,2))
## plot 1st RW
plot.ts(xx, ylab=expression(italic(x[t])))
## plot 2nd RW
plot.ts(x2, ylab=expression(italic(x[t])))

## AUTOREGRESSIVE (AR) MODELS ##
## AR Models of order p, abbreviated AR(p), are commonly used in time series analyses. In particular, ##
## AR(1) models (and their multivariate extensions) see considerable use in ecology as we will see    ##
## later in the course. Recall from lecture that an AR(p) model is written as:                        ##
## xt = φ1xt−1 +φ2xt−2 +···+φpxt−p +wt,  where {wt} is a white noise sequence with zero mean and      ##
## some variance σ2. For our purposes we usually assume that wt ∼ N(0, q).                            ##
## SIMULATE THE AR(p) PROCESS:
## Use the arima.sim function:
set.seed(456)
## list description for AR(1) model with small coef
AR.sm <- list(order=c(1,0,0), ar=0.1, sd=0.1)
## list description for AR(1) model with large coef
AR.lg <- list(order=c(1,0,0), ar=0.9, sd=0.1)
## simulate AR(1)
AR1.sm <- arima.sim(n=50, model=AR.sm)
AR1.lg <- arima.sim(n=50, model=AR.lg)
## Plot the 2 simulated series:
## setup plot region
par(mfrow=c(1,2))
## get y-limits for common plots
ylm <- c(min(AR1.sm,AR1.lg), max(AR1.sm,AR1.lg))
## plot the ts
plot.ts(AR1.sm, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi," = 0.1")))
plot.ts(AR1.lg, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi," = 0.9")))

## Generate two AR(1) models that have the same magnitude coeficient, but opposite signs, and ##
## compare their behavior:
set.seed(123)
## list description for AR(1) model with small coef
AR.pos <- list(order=c(1,0,0), ar=0.5, sd=0.1)
## list description for AR(1) model with large coef
AR.neg <- list(order=c(1,0,0), ar=-0.5, sd=0.1)
## simulate AR(1)
AR1.pos <- arima.sim(n=50, model=AR.pos)
AR1.neg <- arima.sim(n=50, model=AR.neg)
## setup plot region
par(mfrow=c(1,2))
## get y-limits for common plots
ylm <- c(min(AR1.pos,AR1.neg), max(AR1.pos,AR1.neg))
## plot the ts
plot.ts(AR1.pos, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi[1]," = 0.5")))
plot.ts(AR1.neg,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi[1]," = -0.5")))

## Simulate four stationary AR(p) models of increasing order p and then examine their ACF’s and PACF’s:
set.seed(123)
## the 4 AR coefficients
ARp <- c(0.7, 0.2, -0.1, -0.3)
## empty list for storing models
AR.mods <- list()
## loop over orders of p
for(p in 1:4) {
  ## assume SD=1, so not specified
  AR.mods[[p]] <- arima.sim(n=10000, list(ar=ARp[1:p]))
}
## set up plot region
par(mfrow=c(4,3))
## loop over orders of p
for(p in 1:4) {
  plot.ts(AR.mods[[p]][1:50],
          ylab=paste("AR(",p,")",sep=""))
  acf(AR.mods[[p]], lag.max=12)
  pacf(AR.mods[[p]], lag.max=12, ylab="PACF")
}
## NOTE: the ACF for an AR(p) process tails off toward zero very slowly, but the PACF goes to zero for ##
## lags > p. This is an important diagnostic tool when trying to identify the order of p in ARMA(p,q)  ##
## models.

## MOVING-AVERAGE (MA) MODELS ##
## A moving-averge process of order q, or MA(q), is a weighted sum of the current random error plus the ##
## q most recent errors, and can be written as: xt = wt +θ1wt−1 +θ2wt−2 +···+θqwt−q,                    ##
## where {wt} is a white noise sequence with zero mean and some variance σ2; for our purposes we        ##
## usually assume that wt ∼N(0,q). Of particular note is that because MA processes are finite sums of   ##
## stationary errors, they themselves are stationary.                                                   ##
## SIMULATE MA(q) PROCESS:
## 3 examples:
set.seed(123)
## list description for MA(1) model with small coef
MA.sm <- list(order=c(0,0,1), ma=0.2, sd=0.1)
## list description for MA(1) model with large coef
MA.lg <- list(order=c(0,0,1), ma=0.8, sd=0.1)
## list description for MA(1) model with large coef
MA.neg <- list(order=c(0,0,1), ma=-0.5, sd=0.1)
## simulate MA(1)
MA1.sm <- arima.sim(n=50, model=MA.sm)
MA1.lg <- arima.sim(n=50, model=MA.lg)
MA1.neg <- arima.sim(n=50, model=MA.neg)
## setup plot region
par(mfrow=c(1,3))
## plot the ts
plot.ts(MA1.sm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = 0.2")))
plot.ts(MA1.lg,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = 0.8")))
plot.ts(MA1.neg,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = -0.5")))

## 4 examples of ACF and PACF for AR(p) models:
set.seed(123)
## the 4 MA coefficients
MAq <- c(0.7, 0.2, -0.1, -0.3)
## empty list for storing models
MA.mods <- list()
## loop over orders of q
for(q in 1:4) {
  ## assume SD=1, so not specified
  MA.mods[[q]] <- arima.sim(n=1000, list(ma=MAq[1:q]))
}
## Plot the time series, ACF's, and PACF's for the 4 MA(q) models:
## set up plot region
par(mfrow=c(4,3))
## loop over orders of q
for(q in 1:4) {
  plot.ts(MA.mods[[q]][1:50],
          ylab=paste("MA(",q,")",sep=""))
  acf(MA.mods[[q]], lag.max=12)
  pacf(MA.mods[[q]], lag.max=12, ylab="PACF")
}
## NOTE: very little qualitative difference in the realizations of the four MA(q) processes.  ##
## As we saw in lecture and is evident from our examples here, however, the ACF for an MA(q)  ##
## process goes to zero for lags > q, but the PACF tails off toward zero very slowly. This is ## 
## an important diagnostic tool when trying to identify the order of q in ARMA(p,q) models.   ##

## AUTO-REGRESSIVE MOVING-AVERAGE (ARMA) MODELS ##
## Write an ARMA(p,q) as a mixture of AR(p) and MA(q) models, such that:  ##
## xt = φ1xt−1 +φ2xt−2 +···+φpxt−p +wt +θwt−1 +θ2wt−2 +···+θqxt−q,        ##
## and the wt are white noise.                                            ##
## SIMULATE ARMA(p,q) MODELS WITH arima:
set.seed(123)
## ARMA(2,2) description for arim.sim()
ARMA22 <- list(order=c(2,0,2), ar=c(-0.7,0.2), ma=c(0.7,0.2))
## mean of process
mu <- 5
## simulated process (+ mean)
ARMA.sim <- arima.sim(n=10000, model=ARMA22) + mu
## estimate parameters
arima(x=ARMA.sim, order=c(2,0,2))

## Figure out the orders of p and q:
## Search over several possible model forms and see which of them provides the most parsimonious fit to the data.
## Script to loop over possible dimensions of p and q:
## empty list to store model fits
ARMA.res <- list()
## set counter
cc <- 1
## loop over AR
for(p in 0:3) {
  ## loop over MA
  for(q in 0:3) {
    ARMA.res[[cc]] <- arima(x=ARMA.sim,order=c(p,0,q))
    cc <- cc + 1 }
}
## get AIC values for model evaluation
ARMA.AIC <- sapply(ARMA.res,function(x) x$aic)
## model with lowest AIC is the best
ARMA.res[[which(ARMA.AIC==min(ARMA.AIC))]]

## Alternative method to figure out orders of p and q:
## Use auto.arima function in the package forecast: conduct an automatic search over all possible orders of ##
## ARIMA models that you specify:
## find best ARMA(p,q) model
auto.arima(ARMA.sim, start.p=0, max.p=3, start.q=0, max.q=3)