# BNGM


## Aim

Fit neural ordinary differential equations (NODE) by Bayesian neural gradient matching (BNGM) to nonparametrically infer ecological interactions between hare and lynx in the hare-lynx time series.


## Method


## Running BNGM

All modules necessary to perform the analysis can be executed by running the main script:

``` R

## goal: train observation model
system.time(expr={
source("m3_trainModel_o.r");
})
source("m4_plotModel_o.r");

## goal: cross validation on process model
system.time(expr={
source("m9_crossVal_p.r");
})
source("m10_plotCrossVal_p.r");

## goal: train process model 
system.time(expr={
source("m7_trainModel_p.r");
})
source("m8_plotModel_p.r");

```
