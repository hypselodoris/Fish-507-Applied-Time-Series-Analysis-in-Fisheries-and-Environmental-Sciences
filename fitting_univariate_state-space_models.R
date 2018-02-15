# Fish 507: Applied Time-Series Analysis in Fisheries and Environmental Sciences
# Week 3 Problems

# Load utility functions
source("/Volumes/workspaces/nate/public/coding_work/R/functions.R")
# Load packages
packages <- c("MARSS", "stats", "forecast", "datasets")
ipak(packages)

## MARSS package fits multivariate auto-regressive models of the form found in class document corresponding to this lab 
## in equation 1.1:
## xt = Bxt−1 +u+wt where wt ∼ N(0,Q)
## yt = Zxt +a+vt where vt ∼ N(0,R)
## x0 = μ
## NOTE: this equation is in MATRIX FORM!!!
## Example: fit a univariate AR-1 model observed with error:
# xt = bxt−1 +wt where wt ∼ N(0,q)
# yt = xt +vt where vt ∼ N(0,r) 
# x0 = μ
## Model list for this equation:
mod.list=list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )
## Simulate some AR-1 plus error data:
q=0.1; r=0.1; n=100
y=cumsum(rnorm(n,0,sqrt(q)))+rnorm(n,0,sqrt(r))
## Fit with MARSS() using mod.list above:
fit=MARSS(y, model=mod.list)

## Load Nile River flow data
dat=as.vector(Nile)

## 1) Flat level model (i.e. simple average flow with variability around some level μ):
# yt = μ+vt where vt ∼ N(0,r)
## where yt is the river flow volume at year t.
## Write this model as a univariate state-space model:
# xt = 1×xt−1 +0+wt where wt ∼ N(0,0) 
# yt = 1×xt +0+vt where vt ∼ N(0,r) 
# x0 = μ
## Specify model list accordingly:
mod.nile.0 = list(
  B=matrix(1), U=matrix(0), Q=matrix(0),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )
## Fit model with MARSS():
kem.0 = MARSS(dat, model=mod.nile.0)

## 2) Linear trend in flow model
## Looking at the data, we might expect that a declining average river flow would be better.  ##
## In MARSS form, that model would be:
# xt = 1×xt−1 +u+wt where wt ∼ N(0,0)
# yt = 1×xt +0+vt where vt ∼ N(0,r) 
# x0 = μ
## where u is now the average per-year decline in river flow volume. The model
## is specified as follows:
mod.nile.1 = list(
  B=matrix(1), U=matrix("u"), Q=matrix(0),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )
## Fit model with MARSS():
kem.1 = MARSS(dat, model=mod.nile.1)

## 3) Stochastic level model
## Looking at the flow levels, we might suspect that a model that allows the average flow to  ##
## change would model the data better and we might suspect that there have been sudden, and   ##
## anomalous, changes in the river flow level. We will now model the average river flow at    ##
## year t as a random walk, specifically an autoregressive process which means that average   ##
## river flow is year t is a function of average river flow in year t − 1.
# xt = xt−1 +wt where wt ∼ N(0,q)
# yt = xt +vt where vt ∼ N(0,r) 
# x0 = μ
## yt is the river flow volume at year t. xt is the mean level.
## The model is specified as:
mod.nile.2 = list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )
## Fit model with MARSS():
kem.2 = MARSS(dat, model=mod.nile.2)

## 4) Stochastic level model with drift
## Add a drift to term to our random walk; the u in the process model (x) is the drift term.  ##
## This causes the random walk to tend to trend up or down.
# xt = xt−1 +u+wt where wt ∼ N(0,q)
# yt = xt +vt where vt ∼ N(0,r) 
# x0 = μ
## The model is then specified by changing U to indicate that a u is estimated:
mod.nile.3 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0)
## Fit model with MARSS():
kem.3 = MARSS(dat, model=mod.nile.3)

## NOTE: StructTS function also fits the stochastic level model:
fit.sts = StructTS(dat, type="level")
fit.sts
## The estimates from StructTS() will be different (though similar) from MARSS()  ##
## because StructTS() uses x1 = y1, that is the hidden state at t = 1 is fixed to ##
## be the data at t = 1. That is fine if you have a long data set, but would be   ##
## disastrous for the short data sets typical in fisheries and ecology.
## StructTS() is much, much faster for long time series. The example in ?StructTS ##
## is pretty much instantaneous with StructTS() but takes minutes with the EM     ##
## algorithm that is the default in MARSS(). With the BFGS algo- rithm, it is     ##
## much closer to StructTS():
trees <- window(treering, start = 0)
fitts = StructTS(trees, type = "level")
fitem = MARSS(as.vector(trees),mod.nile.2)
fitbf = MARSS(as.vector(trees),mod.nile.2, method="BFGS")

## COMPARE MODELS WITH AIC AND MODEL WEIGHTS
## To get the AIC or AICc values for a model fit from a MARSS fit, use fit$AIC or ##
## fit$AICc. The log-likelihood is in fit$logLik and the number of estimated      ##
## parameters in fit$num.params. For fits from other functions, try AIC(fit) or   ##
## look at the function documentation.                                            ##
## Put the AICc values Nile models together:
nile.aic = c(kem.0$AICc, kem.1$AICc, kem.2$AICc, kem.3$AICc)
## Calculate the AICc minus the minus AICc in our model set and compute the model ##
## weights. ∆AIC is the AIC values minus the minimum AIC value in your model set. ##
delAIC= nile.aic-min(nile.aic)
relLik=exp(-0.5*delAIC)
aicweight=relLik/sum(relLik)
## Create a model weights table:
aic.table=data.frame(
  AICc=nile.aic,
  delAIC=delAIC,
  relLik=relLik,
  weight=aicweight)
rownames(aic.table)=c("flat level","linear trend", "stoc level", "stoc level w drift")
## Print the table with digits limited to specified amount using round():
round(aic.table, digits = 3)
## NOTE: When comparing models within a set of models, at least one of the models must fit the data "reasonably well". ##
## Otherwise, you will be choosing the "less-bad" model among bad models.                                              ##

## BASIC DIAGNOSTICS
## 1) With ANY statistical analysis: compare residuals to assumed error structure...they should correspond. ##
## We have two types of errors in a univariate state-space model: process errors, the wt , and observation  ##
## errors, the vt. They should not have a temporal trend. To get the residuals from most types of fits in R,##
## you can use residuals(fit). MARSS calls the vt, model residuals, and the wt state residuals. We can plot ##
## these using the following code:
par(mfrow=c(3,2))
resids=residuals(kem.0)
plot(resids$model.residuals[1,],
     ylab="model residual", xlab="", main="flat level")
abline(h=0)
plot(resids$state.residuals[1,],
     ylab="state residual", xlab="", main="flat level")
abline(h=0)
resids=residuals(kem.1)
plot(resids$model.residuals[1,],
     ylab="model residual", xlab="", main="linear trend")
abline(h=0)
plot(resids$state.residuals[1,],
     ylab="state residual", xlab="", main="linear trend")
abline(h=0)
resids=residuals(kem.2)
plot(resids$model.residuals[1,],
     ylab="model residual", xlab="", main="stoc level")
abline(h=0)
plot(resids$state.residuals[1,],
     ylab="state residual", xlab="", main="stoc level")
abline(h=0)
## Check the autocorrelation of residuals. THEY SHOULD NOT BE AUTOCORRELATED IN TIME. ##
## No need to check autocorrelation for flat level or linear trends because for those ##
## models wt =0/                                                                      ##
## Calculate acf's:
par(mfrow=c(2,2))
resids=residuals(kem.0)
acf(resids$model.residuals[1,], main="flat level v(t)")
resids=residuals(kem.1)
acf(resids$model.residuals[1,], main="linear trend v(t)")
resids=residuals(kem.2)
acf(resids$model.residuals[1,], main="stoc level v(t)")
acf(resids$state.residuals[1,], main="stoc level w(t)")
## The stochastic level model looks the best in that its model residuals (the vt) are ##
## fine but the state model still has problems. Clearly the state is not a simple     ##
## random walk. This is not surprising. The Aswan Low Dam was completed in 1902 and   ##
## changed the mean flow. The Aswan High Dam was completed in 1970 and also affected  ##
## the flow.

## FITTING A UNIVARIATE AR-1 STATE-SPACE MODEL WITH JAGS*
## Write the model for JAGS to a file (filename in model.loc):
model.loc="ss_model.txt"
jagsscript = cat("
   model {
   # priors on parameters
   mu ~ dnorm(Y1, 1/(Y1*100)); # normal mean = 0, sd = 1/sqrt(0.01)
   tau.q ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.q <- 1/sqrt(tau.q); # sd is treated as derived parameter
   tau.r ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.r <- 1/sqrt(tau.r); # sd is treated as derived parameter
   u ~ dnorm(0, 0.01);
   # Because init X is specified at t=0
   X0 <- mu
   X[1] ~ dnorm(X0+u,tau.q);
   Y[1] ~ dnorm(X[1], tau.r);
   for(i in 2:N) {
   predX[i] <- X[i-1]+u;
   X[i] ~ dnorm(predX[i],tau.q); # Process variation
   Y[i] ~ dnorm(X[i], tau.r); # Observation variation
   }
   }
   ",file=model.loc)
## Specify data and any other input that the JAGS code needs:
## We need:
# dat
# number of time steps 
# parameters we want to monitor
## Note, that the hidden state is a parameter in the Bayesian context (but not in the maximum likelihood context).  ##
jags.data = list("Y"=dat, "N"=length(dat), Y1=dat[1])
jags.params=c("sd.q","sd.r","X","mu", "u")
## Fit the model:
mod_ss = jags(jags.data, parameters.to.save=jags.params,
  model.file=model.loc, n.chains = 3,
  n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)
## Show the posteriors with the MLEs from MARSS on top:
attach.jags(mod_ss)
par(mfrow=c(2,2))
hist(mu)
abline(v=coef(kem.3)$x0, col="red")
hist(u)
abline(v=coef(kem.3)$U, col="red")
hist(log(sd.q^2))
abline(v=log(coef(kem.3)$Q), col="red")
hist(log(sd.r^2))
abline(v=log(coef(kem.3)$R), col="red")
detach.jags()
## Plot the estimated states:
## Custom model output plot function:
plotModelOutput = function(jagsmodel, Y) {
  attach.jags(jagsmodel)
  x = seq(1,length(Y))
  XPred = cbind(apply(X,2,quantile,0.025), apply(X,2,mean), apply(X,2,quantile,0.975)) 
  ylims = c(min(c(Y,XPred), na.rm=TRUE), max(c(Y,XPred), na.rm=TRUE))
  plot(Y, col="white",ylim=ylims, xlab="",ylab="State predictions")
  polygon(c(x,rev(x)), c(XPred[,1], rev(XPred[,3])), col="grey70",border=NA)
  lines(XPred[,2])
  points(Y) 
}

plotModelOutput(mod_ss, dat)
lines(kem.3$states[1,], col="red")
lines(1.96*kem.3$states.se[1,]+kem.3$states[1,], col="red", lty=2)
lines(-1.96*kem.3$states.se[1,]+kem.3$states[1,], col="red", lty=2)
title("State estimate and data from\nJAGS (black) versus MARSS (red)")
## NOTE: Bayesian fit along with 95% credible inter- vals (black and grey)  ##
## MLE states and 95% condidence intervals in red.                          ##

## A SIMPLE RANDOM WALK MODEL OF ANIMAL MOVEMENT
## A random walk model of movement with drift (directional movement) but no correlation:  ##
# x1,t = x1,t−1 +u1 +w1,t, w1,t ∼ N(0,σ21)
# x2,t = x2,t−1 +u2 +w2,t, w2,t ∼ N(0,σ2)
## Where x1,t is the location at time t along one axis (here, longitude) and x2,t is for  ##
## another, generally orthogonal, axis (in here, latitude). The parameter u1 is the rate  ##
## of longitudinal movement and u2 is the rate of latitudinal movement.                   ##
## Add errors to observations of location:
# y1,t = x1,t + v1,t , v1,t ∼ N(0, η21 ) 
# y2,t = x2,t + v2,t , v2,t ∼ N(0, η2 )
## This model is comprised of two separate univariate state-space models. Note that y1    ##
## depends only on x1 and y2 depends only on x2. There are no actual interactions between ##
## these two univariate models. However, we can write the model down in the form of a     ##
## multivariate model using diagonal variance-covariance matrices and a diagonal design   ##
## (Z) matrix. Because the variance-covariance matrices and Z are diagonal, the x1:y1 and ##
## x2:y2 processes will be independent as intended.                                       ##
## Equations written as a MARSS model (in matrix form): See pdf of class materials        ##

## PROBLEMS ##
## 1.1 Write the equations fo reach of these models:ARIMA(0,0,0),ARIMA(0,1,0),            ##
## ARIMA(1,0,0), ARIMA(0,0,1), ARIMA(1,0,1). Read the help file for Arima (forecast       ##
## package) if you are fuzzy on the arima notation.

## 1.2 The MARSS package includes a data set of sharp-tailed grouse in Wash- ington.      ##
## Load the data to use as follows:                                                       ##
library(MARSS)
dat=log(grouse[,2])
## Consider these two models for the data:                                                ##
## Model 1 random walk with no drift observed with no error                               ##
## Model 2 random walk with drift observed with no error                                  ##
## Written as a univariate state-space model, model 1 is                                  ##
# xt = xt−1 +wt where wt ∼ N(0,q) x0 = a
# yt = xt
## Model 2 is almost identical except with u added:
# xt = xt−1 +u+wt where wt ∼ N(0,q) x0 = a
# yt = xt
## y is the log grouse count in year t.
## a) Plot the data. The year is in column 1 of grouse.
plot(dat,
     ylab="Log-Count", xlab="Year", main="Log-Grouse Data")
## b) Fit each model using MARSS().
# The MARSS package fits multivariate auto-regressive models of this form:
#   xt = Bxt−1 +u+wt where wt ∼ N(0,Q)
#   yt = Zxt +a+vt where vt ∼ N(0,R) 
#   x0 = μ
# NOTE: MATRIX FORM! To use MARSS() you must convert a state-space model into the matrix form required by MARSS().
# Define following matrix terms:
# B,U,Q,Z,A,R,x0
# Set which time step x is initialized: tinitx
# Important: in a state-space model, y is always the data and x is something estimated from the data.

## Model 1 as univariate state-space model:
# xt = xt−1 +wt where wt ∼ N(0,q) 
# x0 = a
# yt = xt
## The model list for Model 1:
mod.grouse.1 = list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("a"), tinitx=0)
## Fit model with MARSS():
grouse.marss.fit.1 = MARSS(dat, model = mod.grouse.1)

## Model 2 as univariate state-space model:
# xt = xt−1 +u+wt where wt ∼ N(0,q) 
# x0 = a
# yt = xt
mod.grouse.2 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("a"), tinitx=0)
## Fit model with MARSS():
grouse.marss.fit.2 = MARSS(dat, model = mod.grouse.2)

## COMPARE MODELS WITH AIC AND MODEL WEIGHTS
## To get the AIC or AICc values for a model fit from a MARSS fit, use fit$AIC or ##
## fit$AICc. The log-likelihood is in fit$logLik and the number of estimated      ##
## parameters in fit$num.params. For fits from other functions, try AIC(fit) or   ##
## look at the function documentation.                                            ##
## Put the AICc values Nile models together:
grouse.marss.aic = c(grouse.marss.fit.1$AICc, grouse.marss.fit.2$AICc)
## Calculate the AICc minus the minus AICc in our model set and compute the model ##
## weights. ∆AIC is the AIC values minus the minimum AIC value in your model set. ##
delAIC= grouse.marss.aic-min(grouse.marss.aic)
relLik=exp(-0.5*delAIC)
aicweight=relLik/sum(relLik)
## Create a model weights table:
aic.marss.table=data.frame(
  AICc=grouse.marss.aic,
  delAIC=delAIC,
  relLik=relLik,
  weight=aicweight)
rownames(aic.marss.table)=c("model 1","model 2")
## Print the table with digits limited to specified amount using round():
round(aic.marss.table, digits = 3)
## Output:
#           AICc delAIC relLik weight
# model 1  0.634  2.568  0.277  0.217
# model 2 -1.934  0.000  1.000  0.783
## c) Which one appears better supported given AICc?
# Answer: "model 2" appears better supported.

## d) Load the forecast package. Use ?auto.arima to learn what it does. ##
## Then use auto.arima(dat) to fit the data. Next run auto.arima on the ##
## data with trace=TRUE to see all the ARIMA models it compared. Note,  ##
## ARIMA(0,1,0) is a random walk with b=1. ARIMA(0,1,0) with drift      ##
## would be a random walk (b=1) with drift (with u).
## OUTPUT:
# auto.arima(dat, trace = TRUE)
# 
# ARIMA(2,1,2) with drift         : Inf
# ARIMA(0,1,0) with drift         : -3.116841
# ARIMA(1,1,0) with drift         : -1.006099
# ARIMA(0,1,1) with drift         : Inf
# ARIMA(0,1,0)                    : -0.5520726
# ARIMA(1,1,1) with drift         : Inf
# 
# Best model: ARIMA(0,1,0) with drift         
# 
# Series: dat 
# ARIMA(0,1,0) with drift 
# 
# Coefficients:
#   drift
# -0.0909
# s.e.   0.0394
# 
# sigma^2 estimated as 0.0467:  log likelihood=3.79
# AIC=-3.58   AICc=-3.12   BIC=-0.84
## NOTE: It picked model 2 as the best among those tested. ”ARIMA(0,1,0) with drift” is model 2.

## e) Is the difference in the AICc values between a random walk with ##
## and without drift comparable between MARSS() and auto.arima()?
grouse.arima.fit.1.arima = Arima(dat, order=c(0,1,0))
grouse.arima.fit.2.arima = Arima(dat, order=c(0,1,0), include.drift=TRUE)
grouse.arima.fit.2.arima$aicc - grouse.arima.fit.1.arima$aicc
## OUTPUT:
# [1] -2.564768
grouse.marss.fit.2$AICc - grouse.marss.fit.1$AICc
## OUTPUT:
# [1] -2.567739
## NOTE: Similar, but not identical results.
## When using auto.arima, it will fit AR-1 models of the following form (notice the b): ##
## xt = bxt−1 +wt. auto.arima refers this model xt = xt−1 +wt, which is also AR-1 but   ##
## with b = 1, as ARIMA(0,1,0). This says that the first difference of the data (that’s ##
## the 1 in the middle) is a ARMA(0,0) process (the 0s in the 1st and 3rd spots). So    ##
## ARIMA(0,1,0) means this: xt −xt−1 =wt.

## 3) Create a random walk with drift time series using the following command:          ##
## dat=cumsum(rnorm(100,0.1,1))
## NOTE: n = 100, mean = 0.1, sd = 1
## a) Write out the equation for that random walk as a univariate state-space model.    ##
# xt = xt−1 +u+wt,wt ∼ N(0,q)
# x0 =μ or x1 =y1 (1.1)
# yt = xt
# where u=0.1 and q=1.

## b) What is the order of the x part of the model written as ARIMA(p, d, q)?           ##
# From question 1, you should be able to deduce it is ARIMA(0,1,0) but if you said      ##
# ARIMA(1,0,0) with b=1, that’s ok. That’s not how Arima() writesxt=xt−1+u+wt but it is ##
# correct.

## c) Fit that model using Arima() in the forecast package. You’ll need to specify the  ##
## order and include.drift term.                                                        ##
arima.fit.3.d = Arima(dat, order=c(0,1,0), include.drift=TRUE)

## d) Fit that model with MARSS().                                                      ##
mod.3.d = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("μ"), tinitx=0)
## Fit model with MARSS():
marss.fit.3.d = MARSS(dat, model = mod.3.d)

## e) How are the two estimates different?
coef(marss.fit.3.d, type="vector")
c(coef(arima.fit.3.d), s2=arima.fit.3.d$sigma2)
## Now fit the first-differenced data:
diff.dat=diff(dat)

## f) If xt denotes a time series. What is the first difference of x? What is the     ##
## second difference?                                                                 ##
# First difference diff(x) is xt −xt−1.
# Second difference is diff(diff(x)) or (xt −xt−1)−(xt−1 −xt−2).

## g) What is the x model for diff.dat?
# diff(x) = (xt −xt−1) = u+wt

## h) Fit diff.dat using Arima().You’ll need to change order and include.mean.        ##
fit.diff.Arima=Arima(diff.dat, order=c(0,0,0), include.mean=TRUE)
fit.diff.arima=arima(diff.dat, order=c(0,0,0), include.mean=TRUE)

## i) Fit that model with MARSS()
# data (y) is now diff.dat and state-space model is
# xt = u+wt,wt ∼ N(0,q) x0 = 0
# yt = xt
## It doesn’t matter what x0 is; it does not appear in the model, but it is important to use x0      ##
## instead of x1 to match arima().                                                                   ##
mod.diff.dat = list(
  B=matrix(0), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix(0), tinitx=0)
fit.diff.marss = MARSS(diff.dat, model = mod.diff.dat)

## Note, we can also fit with lm():                                                                 ##                                                  
fit.diff.lm = lm(diff.dat~1)
## Here are the parameter estimates:                                                                ## 
rbind(
  marss.diff=coef(fit.diff.marss, type="vector"), 
  arima.diff=c(coef(fit.diff.arima), s2=fit.diff.arima$sigma2), 
  Arima.diff=c(coef(fit.diff.Arima), s2=(98/99)*fit.diff.Arima$sigma2), 
  lm.diff=c(coef(fit.diff.lm), s2=(98/99)*summary(fit.diff.lm)$sigma^2)
)
## OUTPUT:                                                                                          ##
#                   U.u      Q.q
# marss.diff 0.03824377 1.001234
# arima.diff 0.03824377 1.001234
# Arima.diff 0.03824377 1.001234
# lm.diff    0.03824377 1.001234
## They are all the same except the variances reported by Arima() and lm() have to be multiplied by ##
## 98/99 to be the same as MARSS and arima because the former are reporting the unbiased estimates  ##
## and the latter are reporting the straight (biased) maximum-likelihood estimates.                 ##

## 4) Arima() will also fit what it calls an ‘AR-1 with drift’. An AR-1 with drift is NOT this model: ##
# xt = bxt−1 +u+wt where wt ∼ N(0,q)    <- Equation 1.4
## In the population dynamics literature, this equation is called the Gompertz model and is a type  ##
## of density-dependent population model.                                                           ##
## a) Write R code to simulate Equation 1.4. Make b less than 1 and greater                         ##
## than 0. Set u and x0 to whatever you want. You can use a for loop.                               ##
# set up my parameter values:
b=.8; u=2; x0=10; q=0.1
nsim=1000
# set up my holder for x:
x=rep(NA, nsim)
# first time step:
x[1]=b*x0+u+rnorm(1,0,sqrt(q))
# use for loop for all time steps after t=1
for(t in 2:nsim) x[t]=b*x[t-1]+u+rnorm(1,0,sqrt(q))

## b) Plot the trajectories and show that this model does not “drift” upward or downward. It        ##
## fluctuates about a mean value.                                                                   ##
plot(x, type="l",xlab="", ylab="x")

## c) Hold b constant and change u. How do the trajectories change?                                 ##
# Change u, add 1:
u2=u+1
# Set up new holder for x2:
x2=rep(NA, nsim)
# First time step for x2:
x2[1]=b*x0+u2+rnorm(1,0,sqrt(q))
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x2[t]=b*x2[t-1]+u2+rnorm(1,0,sqrt(q))

# Change u, subtract 1:
u3=u-1
# Set up new holder for x3:
x3=rep(NA, nsim)
# First time step for x3:
x3[1]=b*x0+u3+rnorm(1,0,sqrt(q))
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x3[t]=b*x3[t-1]+u3+rnorm(1,0,sqrt(q))

# Set up plot device: 3 rows, 1 column
par(mfrow=c(3,1))
plot(x, type="l", main="X", xlab="", ylab="x")
plot(x2, type="l", main="X2", xlab="", ylab="x")
plot(x3, type="l", main="X3", xlab="", ylab="x")
# Answer: Means are shifted up and down respectively. 

## d) Hold u constant and change b. Make sure to use a b close to 1 and   ##
## another close to 0. How do the trajectories change?                    ##
# set up my parameter values:
# Set b2 close to 1:
b2=.9
# Set up new holder for x4:
x4=rep(NA, nsim)
# First time step for x4:
x4[1]=b2*x0+u+rnorm(1,0,sqrt(q))
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x4[t]=b2*x4[t-1]+u+rnorm(1,0,sqrt(q))

# Set b3 close to 0:
b3=.1
# Set up new holder for x5:
x5=rep(NA, nsim)
# First time step for x5:
x5[1]=b3*x0+u+rnorm(1,0,sqrt(q))
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x5[t]=b3*x5[t-1]+u+rnorm(1,0,sqrt(q))

# Plot trajectories of original and x4, x5:
range(x)
range(x4)
range(x5)

# plot(x, type="l", ylim=range(1,23), xlab="", ylab="x")
plot(x4, type="l", ylim=range(1,23), col="black")
lines(x5, col="red")
legend('topright', legend = c("b=0.9", "b=0.1"), col = c("black", "red"), lty =1)
# The one with smaller b has less auto-regression and is ‘tighter’ (explores less of a range of the y axis).

## e) Do 2 simulations each with the same wt . In one simulation, set u = 1 and in the other u = 2. ##
## For both simulations, set x1 = u/(1−b). You can set b to whatever you want as long as 0 < b < 1. ## 
## Plot the 2 trajectories on the same plot. What is different?                                     ##
# Trajectory 1:
nsim = 1000
b = 0.8
u = 1
x0 = u/(1-b)
q = 0.1
err = rnorm(nsim, 0, q)
# Set up variable to hold simulated data:
x = rep(NA, nsim)
# First time step for x:
x[1] = b * x0 + u + err[1]
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x[t] = b * x[t-1] + u + err[t]

# Trajectory 2:
u = 2
# Set up variable to hold simulated data:
x2 = rep(NA, nsim)
# First time step for x:
x2[1]=b*x0+u+err[1]
# Use for loop to iterate over all time steps after t=1
for(t in 2:nsim) x2[t]=b*x2[t-1]+u+err[t]

# Plot both simulations:
plot(x, type="l", ylim=c(3.5, 12.5), col="black")
lines(x2, col="red")
legend("topright", legend = c("u=1", "u=2"), col = c("black", "red"), lty = 1)
## NOTE: Trajectories 1 and 2 are identical except that the mean of trajectory 2 is shifted up.                 ##
## NOTE: We will fit what Arima calls “AR-1 with drift” models in the chapter on MARSS models with covariates.  ##

## 5) The MARSS package includes a data set of gray whales. Load the data to use as follows:
library(MARSS)
dat=log(graywhales[,2])
## Fit a random walk with drift model observed with error to the data:
# xt = xt−1 +u+wt where wt ∼ N(0,q)
# yt = xt +vt where vt ∼ N(0,r) 
# x0 = a
## y is the whale count in year t. x is interpreted as the ’true’ unknown population size that we are trying  ##
## to estimate.                                                                                               ##
## a) Fit this model with MARSS()                                                                             ##
mod.whale = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("a"), tinitx=0
)
whale.marss.fit = MARSS(dat, model = mod.whale)
## OUTPUT:
# Success! abstol and log-log tests passed at 16 iterations.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Estimation converged in 16 iterations. 
# Log-likelihood: 4.064946 
# AIC: -0.129891   AICc: 1.975372   
# 
# Estimate
# R.r    0.0141
# U.u    0.0564
# Q.q    0.0136
# x0.a   7.9532
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.

## b) Plot the estimated x as a line with the actual counts added as points.x is in whale.marss.fit$states.   ##
## It is a matrix, which plot will not like so you will need to change it to a vector using as.vector() or    ##
## whale.marss.fit$states[1,]
x = as.vector(whale.marss.fit$states)
plot(dat, type = "p", main = "Log-Gray Whales", ylab = "Log-counts", xlab = "", col = "black")
lines(x, col = "red")
# Plot with years labelled:
par(mar=c(2,2,2,2))
plot(graywhales[,1], whale.marss.fit$states[1,], type="l",xlab="", ylab="log count") 
points(graywhales[,1], dat)


## c) Simulate 1000 sample trajectories using the estimated u and q starting at the estimated x in 1997. You  ##
## can do this with a couple for loops or write something terse with cumsum and apply.                        ##
n.trajectories = 1000
n.forward = 10
b = 0.8
u = as.vector(coef(whale.marss.fit)$U)
x0 = whale.marss.fit$states[1,39]
q = as.vector(coef(whale.marss.fit)$Q)
# Set up matrix to hold simulated data. Each row contains 1 trajectory, each column is a time step foreward:
x = matrix(NA, n.trajectories, n.forward)
# Set up variable to hold simulated data for all trajectories:
# First time step for x:
x[,1]=x0+u+rnorm(n.trajectories,0,sqrt(q))
# Iterate over number of years foreward, populate data matrix:
for(t in 2:n.forward) x[,t]=x[,t-1]+u+rnorm(n.trajectories,0,sqrt(q))

## d) Using these simulated trajectories, what is your estimated probability of reaching 50,000 graywhales  ##
## in 2007.
# I just want the fraction of simulations that were 50,000 or above in 2007
x.threshold = log(50000)
sum(x[,10]<=x.threshold)/n.trajectories

## e) What kind of uncertainty does that estimate NOT include?                                              ##
# By using the point estimates of u, q and x0, we are not including the uncertainty in those estimates in   ##
# our forecasts.

## 6) Fit the following models to the graywhales data using MARSS(). Assume b=1.                            ##
## Model 1 Process error only model with drift:
# xt = xt−1 +u+wt where wt ∼ N(0,q)
# yt = xt
# x0 = a
mod.list.1 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0
)
fit.marss.1 = MARSS(dat, model = mod.list.1)

## Model 2 Process error only model without drift:
# xt = xt−1 +wt where wt ∼ N(0,q)
# yt = xt
# x0 = mu
mod.list.2 = list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0
)
fit.marss.2 = MARSS(dat, model = mod.list.2)

## Model 3 Process error with drift and observation error with observation error variance fixed = 0.05:
# xt = xt−1 +u+wt where wt ∼ N(0,q)
# yt = xt +vt where vt ~ N(0,0.05)
# x0 = mu
mod.list.3 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0.05),
  x0=matrix("mu"), tinitx=0
)
fit.marss.3 = MARSS(dat, model = mod.list.3)

## Model 4 Process error with drift and observation error with observation error variance estimated:
# xt = xt−1 +u+wt where wt ∼ N(0,q)
# yt = xt +vt where vt ~ N(0,r)
# x0 = mu
mod.list.4 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0
)
fit.marss.4 = MARSS(dat, model = mod.list.4)

## a) Compute the AICc’s for each model and likelihood or deviance (-2 * log likelihood). Where to find ##
## these? Try names(fit). logLik() is the standard R function to return log-likelihood from fits.       ##
fit.aic = c(fit.marss.1$AICc, fit.marss.2$AICc, fit.marss.3$AICc, fit.marss.4$AICc)
fit.logLik = c(fit.marss.1$logLik, fit.marss.2$logLik, fit.marss.3$logLik, fit.marss.4$logLik)
## Create a table with AICc and logLik:
aic.table=data.frame(
  AICc=fit.aic,
  logLik=fit.logLik
  )
rownames(aic.table)=c("Model 1","Model 2", "Model 3", "Model 4")
## Print the table with digits limited to specified amount using round():
round(aic.table, digits = 3)

## b) Calculate a table of ∆AICc values and AICc weights.                                               ##
# Calculate ∆AICc values and AICc weights:
delAIC= fit.aic-min(fit.aic)
relLik=exp(-0.5*delAIC)
aicweight=relLik/sum(relLik)
# Build a table of results:
aic.table=data.frame(
  AICc=fit.aic,
  logLik=fit.logLik,
  delAICc=delAIC,
  relLik=relLik,
  weight=aicweight
)
rownames(aic.table)=c("Model 1","Model 2", "Model 3", "Model 4")
## Print the table with digits limited to specified amount using round():
round(aic.table, digits = 3)

## c) Show the acf of the model and state residuals for the best model. You will need a vector of the   ##
## residuals to do this. If fit is the fit from a fit call like fit = MARSS(dat), you get the residuals ##
## using this code:
# residuals(fit)$state.residuals[1,]
# residuals(fit)$model.residuals[1,]
## Do the acf’s suggest any problems?                                                                   ##
## NOTE: Check the autocorrelation of residuals. THEY SHOULD NOT BE AUTOCORRELATED IN TIME.             ##
model.resids=residuals(fit.marss.4)$model.residuals[1,]
state.resids=residuals(fit.marss.4)$state.residuals[1,]

par(mfrow=c(2,1))
acf(model.resids, main="model residuals")
acf(state.resids, main="state residuals")
# 1 significant point at lag=1, but not concerning because it could be expected by chance alone.

## 7) Evaluate the predictive accuracy of forecasts using the forecast package using the ’airmiles’     ##
## dataset. Load the data to use as follows:                                                            ##
library(forecast)
dat=log(airmiles)
n=length(dat)
# This will prepare the training data and set aside the last 3 data points for validation.:
training.dat = dat[1:(n-3)]
test.dat = dat[(n-2):n]

## a) Fit the following four models using Arima():ARIMA(0,0,0),ARIMA(1,0,0), ARIMA(0,0,1), ARIMA(1,0,1).
# ARIMA(0,0,0) :
arima.fit.1 = Arima(training.dat, order = c(0,0,0))
# ARIMA(1,0,0) :
arima.fit.2 = Arima(training.dat, order = c(1,0,0))
# ARIMA(0,0,1) :
arima.fit.3 = Arima(training.dat, order = c(0,0,1))
# ARIMA(1,0,1) :
arima.fit.4 = Arima(training.dat, order = c(1,0,1))

## b) Use forecast() to make 3 step ahead forecasts from each.
forecast.1 = forecast(arima.fit.1, h = 3)
forecast.2 = forecast(arima.fit.2, h = 3)
forecast.3 = forecast(arima.fit.3, h = 3)
forecast.4 = forecast(arima.fit.4, h = 3)

## c) Calculate the MASE statistic for each using the accuracy function in the forecast package.  ##
## Type ?accuracy to learn how to use this function.                                              ##
accuracy(forecast.1, test.dat)
# The MASE statistic we want is in the Test set row and MASE column.
# Build a numeric vector with results of interest:
MASE.stats = c(
  accuracy(forecast.1, test.dat)["Test set", "MASE"],
  accuracy(forecast.2, test.dat)["Test set", "MASE"],
  accuracy(forecast.3, test.dat)["Test set", "MASE"],
  accuracy(forecast.4, test.dat)["Test set", "MASE"]
)

## d) Present results in a table:
model.names = c("ARIMA(0,0,0)", "ARIMA(1,0,0)", "ARIMA(0,0,1)", "ARIMA(1,0,1)")
MASE.table=data.frame(
  model.name=model.names,
  MASE=MASE.stats
)
# Render table:
MASE.table

## e) Which model is best supported based on the MASE statistic?
# The ARIMA(1,0,1) has the lowest MASE (Mean Absolute Scaled Error) statistic, and is therefore ##
# the best. NOTE:  the AR component strongly improves predictions                               ##

## 8) The WhaleNet Archive of STOP Data has movement data on loggerhead turtles on the east     ##
## coast of the US from ARGOS tags. The MARSS package loggerheadNoisy dataset is lat/lot data   ##
## on eight individuals, however we have corrupted this data severely by adding random errors   ##
## in order to create a “bad tag” problem (very noisy). Use head(loggerheadNoisy) to get an     ##
## idea of the data. Then load the data on one turtle, MaryLee. MARSS needs time across the     ##
## columns so you need to use transpose the data:
turtlename="MaryLee"
dat = loggerheadNoisy[which(loggerheadNoisy$turtle==turtlename),5:6]
dat = t(dat)

## a) Plot MaryLee’s locations (as a line not dots). Put the latitude lo- cations on the y-axis ##
## and the longitude on the y-axis. You can use rownames(dat) to see which is in which row. You ##
## can just use plot() for the homework. But if you want, you can look at the MARSS Manual      ##
## chapter on animal movement to see how to plot the turtle locations on a map using the maps   ##
## package.
#load the map package; you have to install it first
library(maps)
# Read in our noisy data (no missing values)
pdat = loggerheadNoisy #for plotting
turtlename="MaryLee"
par(mai = c(0,0,0,0),mfrow=c(1,1))
map('state', region = c('florida', 'georgia', 'south carolina', 'north carolina',
  'virginia', 'delaware','new jersey','maryland'),xlim=c(-85,-70))
points(pdat$lon[which(pdat$turtle==turtlename)], pdat$lat[which(pdat$turtle==turtlename)],
  col="blue",pch=21, cex=0.7)
lines(pdat$lon[which(pdat$turtle==turtlename)], pdat$lat[which(pdat$turtle==turtlename)],
  col="blue")

## b) Analyze the data with a state-space model (movement observed with error) using:            ##
fit0 = MARSS(dat)
# OUTPUT:
# Success! algorithm run for 15 iterations. abstol and log-log tests passed.
# Alert: conv.test.slope.tol is 0.5.
# Test with smaller values (<0.1) to ensure convergence.
# 
# MARSS fit is
# Estimation method: kem 
# Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
# Algorithm ran 15 (=minit) iterations and convergence was reached. 
# Log-likelihood: -126.2077 
# AIC: 266.4154   AICc: 267.0901   
# 
# Estimate
# R.diag            0.0948  # Rdiag is the observation error variance.
# U.X.lon           0.0753  # Ulon is the average velocity in N-S direction.
# U.X.lat           0.0726  # Ulat is the average velocity in E-W direction.
# Q.(X.lon,X.lon)   0.1024  # Qlon: movement error variance in N-S direction.
# Q.(X.lat,X.lat)   0.0915  # Qlat: movement error variance in E-W direction.
# x0.X.lon        -81.2401  # x0lon: estimated longitudinal position at t=0.
# x0.X.lat         31.7818  # x0lat: estimated latitudinal position at t=0.
# 
# Standard errors have not been calculated. 
# Use MARSSparamCIs to compute CIs and bias estimates.

## c) What assumption did the default MARSS model make about observation error and process error? ##
# The observation errors in the lat and lon direction are independent but have identical variance.##
# The movement errors are independent (not corre- lated) and allowed to have different variances. ##
# So the model doesn’t allow a average NE movement; that would require correlation in the movement##
# errors. It allows that turtles tend to move faster N-S (along the coast) than E-W (out to sea). ##

## d) Does MaryLee move faster in the latitude direction versus longitude direction?              ##
# No. The estimated u’s in the lat and lon direction are similar.

## e) Add MaryLee’s estimated ”true”positions to your plot of her locations. You can use          ##
## lines(x, y, col="red") (with x and y replaced with your x and y data). The true position is    ##
## the ”state”. This is in the states element of an output from MARSS fit0$states.                ##
lines(fit0$states["X.lon",], fit0$states["X.lat",],
    col="red")

## f) Compare the following models for these data. Movement in the lat/lon direction is (1)       ##
## independent but the variance is the same, (2) is correlated and lat/lon variances are different##
## , and (3) is correlated and the lat/lon variances are the same. You only need to change Q      ##
## specification. Your MARSS call will now look like the following with ... replaced with your Q  ##
## specification: fit1 = MARSS(dat, list(Q=...))
fit1 = MARSS(dat, model = list(Q="diagonal and equal"))
fit2 = MARSS(dat, model = list(Q="unconstrained"))
fit3 = MARSS(dat, model = list(Q="equalvarcov"))

c(fit0$AICc,fit1$AICc, fit2$AICc, fit3$AICc)

## g) Plot your state residuals (true location residuals). What are the problems? Discuss in      ##
## reference to your plot of the location data.                                                   ##
resids.3.lon = residuals(fit3)$state.residuals[1,]
resids.3.lat = residuals(fit3)$state.residuals[2,]
par(mfrow=c(2,2),mar=c(3,5,3,5))
plot(resids.3.lon, xlab="") ; abline(h=0)
acf(resids.3.lon)
plot(resids.3.lat, xlab="") ; abline(h=0)
acf(resids.3.lat)
## NOTE: There is a period in the middle of the track where the model does not describe the       ##
## movement well. We can see in the plot that the turtle has a long northward movement in the     ##
## middle of the track.
