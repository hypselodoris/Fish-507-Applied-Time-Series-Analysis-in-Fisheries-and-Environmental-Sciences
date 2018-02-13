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
