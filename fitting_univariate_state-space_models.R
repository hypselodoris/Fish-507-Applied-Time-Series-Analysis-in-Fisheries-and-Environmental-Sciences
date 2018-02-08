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


