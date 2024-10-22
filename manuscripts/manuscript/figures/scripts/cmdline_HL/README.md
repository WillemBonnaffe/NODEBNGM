# BNGM


## Aim

Fit neural ordinary differential equations (NODE) by Bayesian neural gradient matching (BNGM) to nonparametrically infer ecological interactions between the variables in the time series.

The data time series consists of 90 years of pelt counts of hare and lynx in the Hudson bay area. 
We used the dataset from our previous study (Bonnaffe et al. 2021).

The outcome of the approach is nonparametric estimates of the per-capita growth rate of lynx and hare, as well as mutual effects and contributions to the dynamics of these two species.

Bonnaffé, W., Sheldon, B.C. and Coulson, T., 2021. Neural ordinary differential equations for ecological and evolutionary time‐series analysis. Methods in Ecology and Evolution, 12(7), pp.1301-1315.

### Framework inputs:
* time series of population densities (strictly positive)
* parameters of the observation (interpolation) and process model (NODE)

### Framework outputs:
* interpolated states and dynamics (temporal change in population densities)
* effects and contributions of variable dynamics on each other

## Installation

Requirements:
* R (v4.0.2 or later)


## Initialising BNGM

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

In the script `m6_loadModel_p.r`, specify the parametes of the process model:

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

## goal: cross-validation on process model
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


## Description of the files

### Initialisation files: 
* `f_NODE_GM.r` contains all the functions necessary to running the scripts
* `m0_main` contains the commands to run the scripts in the right order to fit NODEs by BNGM 

### Data:
* `data` contains the time series data as csv files
* `data/TS_LV.csv` contains the time series of the hare and lynx

### Observation model:
* `m1_loadData_o.r` loads and formats the time series data
* `m2_loadModel_o.r` loads the observation model (i.e. neural networks that interpolate the time series) 
* `m3_trainModel_o.r` trains the observation model 
* `m4_plotModel_o.r` visualises outputs of the observation model (i.e. interpolations of states and dynamics) 

### Process model:
* `m5_loadData_p.r` loads and formats the interpolated time series and dynamics 
* `m6_loadModel_p.r` loads the process model (i.e. neural networks that approximates dynamics as function of interpolated states) 
* `m7_trainModel_p.r` trains the process model 
* `m8_plotModel_p.r` visualises outputs of the process model (i.e. model fit, effects, and contributions)
* `m9_crossVal_p.r` performs k-fold cross-validation of the process model 
* `m10_plotCrossVal_p.r` visualises outputs of k-fold cross-validation of the process model

### Output files:
* `out` contains the output of the scripts 
* `out_repeat` contains repeat results of the analysis to assess repeatability of the results 
* `out_threefold` contains repeat results obtained with threefold cross-validation (instead of twofold)
* `out/contribTable.csv` contains the contributions of each variable to the dynamics 
* `out/crossVal_p-Omega_p.RData` contains parameters of the process model for each value of the regularisation hyperparameter 
* `out/crossVal_p.RData` contains the likelihood of the predictions of the model for each value of the regularisation hyperparameter
* `out/ddt.Yhat_o.RData` contains the interpolated dynamics of the state variables 
* `out/ddx.Yhat_p.RData` contains the sensitivity of the per-capita growth rate of variables with respect to the each state variable
* `out/effectsTable.csv` contains the effects of each variable on the dynamics 
* `out/fig_crossVal_p.pdf` displays the cross-validation results 
* `out/fig_predictions_o.pdf` displays results of the observation model, i.e. the interpolation of states and dynamics of each variable
* `out/fig_predictions_p.pdf` displays results of the process model, i.e. the effects and contribution of each variable to the dynamics of the system 
* `out/Geber_p.RData` contains the contributions of each variable to the dynamics of the system 
* `out/Omega_o.RData` contains the ensemble of parameters of the observation model 
* `out/Omega_p.RData` contains the ensemble of parameters of the process model 
* `out/summaryTable.csv` contains the summary statistics of the fit of the process model (NODE) for each variable 
* `out/summaryTable.RData` contains the summary statistics of the fit of the process model (NODE) for each variable 
* `out/Yhat_o.RData` contains the interpolated state variables 
* `out/Yhat_p.RData` contains the interpolated per-capita growth rate of each state variable
