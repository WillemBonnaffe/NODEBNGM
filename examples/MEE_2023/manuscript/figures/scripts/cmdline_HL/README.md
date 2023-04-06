# BNGM


## Aim

Fit neural ordinary differential equations (NODE) by Bayesian neural gradient matching (BNGM) to nonparametrically infer ecological interactions between the variables in the time series.


## Method


## Installation

Requirements:
* R v4.2.0 or later


## Initialising BNGM

The BNGM analysis should be initialised by specifying the inputs to the model.

In the script `m1_loadData_o.r`, specify the path to the time series: 

``` R
TS = read.table("data/TS_HL.csv",sep=";",header=T)
```

In the script `m2_loadModel_o.r`, specify the following parameters of the observation model:

``` R
## parameters observation model
K_o      = 100 # number of ensemble samples
W_o      = rep(100,N) # number of nodes in the hidden layer of the networks 
sd2_o    = rep(0.01,N) # standard deviation of the prior distribution of the network parameters
```

In the script `m2_loadModel_p.r`, specify the parametes of the observation model:

``` R
## parameters process model
train_split = 2/3 # train validation split proportion
K_p         = 10 # number of ensemble samples
W_p         = rep(10,N) # number of nodes in hidden layer of network
sd2_p       = as.list(rep(0.1,N)) # standard deviation of the prior distribution of the network parameters
```


## Running BNGM

When the initialisation is complete, all modules necessary to perform the analysis can be executed by running the main script:

``` R
## goal: train observation model
source("m3_trainModel_o.r")
source("m4_plotModel_o.r")

## goal: cross validation on process model
source("m9_crossVal_p.r")
source("m10_plotCrossVal_p.r")

## goal: train process model 
source("m7_trainModel_p.r")
source("m8_plotModel_p.r")
```

This can be done simply by running the command:
``` bash
Rscript m0_main.r
```
