# Fish 507: Applied Time-Series Analysis in Fisheries and Environmental Sciences
# Fitting multivariate state-space models with covariates
# Week 4 Problems

# Load utility functions
source("/Volumes/workspaces/nate/public/coding_work/R/functions.R")
# Load packages
packages <- c("MARSS", "stats", "forecast", "datasets", "coda", "R2jags", "maps")
ipak(packages)

## A MARSS model with covariate effects in both the process and observation components is written as: ##
# xt = Btxt−1 +ut +Ctct +wt, where wt ∼ MVN(0,Qt)
# yt = Ztxt +at +Dtdt +vt, where vt ∼ MVN(0,Rt)
## where ct is the p×1 vector of covariates (e.g., temperature, rainfall) which affect the states and ##
## dt is a q×1 vector of covariates (potentially the same as ct ), which affect the observations. Ct  ##
## is an m × p matrix of coefficients relating the effects of ct to the m×1 state vector xt, and Dt is##
## an n×q matrix of coefficients relating the effects of dt to the n×1 observation vector yt.         ##

## 1.1 Examples using plankton data                                                                   ##
## Here we show some examples using the Lake Washington plankton data set and covariates in that      ##
## dataset. We use the 10 years of data from 1965-1974 (Figure 1.1), a decade with particularly high  ##
## green and bluegreen algae levels. We use the transformed plankton dataset which has 0s replaced    ##
## with NAs. Below, we set up the data and z-score the data. The original data were already z-scored, ##
## but we changed the mean when we subsampled the years so need to z-score again.
fulldat = lakeWAplanktonTrans
years = fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
dat = t(fulldat[years,c("Greens", "Bluegreens")])
# z-score the states:
the.mean = apply(dat,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(dat,1,var,na.rm=TRUE))
dat = (dat-the.mean)*(1/the.sigma)

## Set up the covariate data, temperature and total phosphorous. We z-score the covariates to         ##
## standardize and remove the mean:
covariates = rbind(
  Temp = fulldat[years,"Temp"],
  TP = fulldat[years,"TP"])
# z.score the covariates:
the.mean = apply(covariates,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates = (covariates-the.mean)*(1/the.sigma)

## 1.2 Observation-error only model                                                                   ##
