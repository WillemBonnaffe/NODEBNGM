# NODEBNGM analysis of Ushio system


## Aim

The aim of this analysis is to estimate the effects and contributions of each species to the dynamics of the populations in the Maizuru bay community (Ushio et al. 2018).
The system consists in 12 years long time series of fortnight abundance estimates of 15 of dominant species of the aquatic community of the Maizuru bay in Japan.
The time series are analysed by fitting neural ordinary differential equations (NODE) via Bayesian neural gradient matching (BNGM), which provide interpolations of the state and dynamics of the species, as well as nonparametric estimates of the per-capita growth rates of the species as a function of species density. 
By computing the sensitivity of the per-capita growth rates with respect to each species density, we were able to derive effects and contributions to the dynamics of the species on each other.

Ushio, M., Hsieh, C.H., Masuda, R., Deyle, E.R., Ye, H., Chang, C.W., Sugihara, G. and Kondoh, M., 2018. Fluctuating interaction network and time-varying stability of a natural fish community. Nature, 554(7692), pp.360-363.


## Inputs

### Data

The first input is the time series data, formatted as a csv file, which contains abundance estimates of the species (counts) and environmental variables (sea-bottom temperature in degrees celsius here).
All variables need to be strictly positive.

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_time_series.png)


### Parameters of the observation process model

The second input is the parameters of the observation model (i.e. the neural networks that interpolates the time series) and process model (the NODEs that approximate the per-capita growth rate based on the interpolated variables).
See below for instructions on how to specify the parameters and run the code.


## Outputs

### Interpolated time series and dynamics

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_predictions_o.png)


### Effects and contributions of variables to the dynamics of other variables

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_predictions_p.png)


### Dynamical interaction netwrosk

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_DIN_v1.png)


## Installation

The approach can be used simply by installing R (v4.0.2 or later).


## Running the script

## Preparing the data

The first

## Fitting process model

### Goal 

Fit process model (i.e. explain the per-capita growth rate of the populations calculated as $1/Y*dY/dt$ as a function of the states Y(t))

### Notes 

* the user could use state interpolations and interpolated dynamics obtained via other methods (e.g. Fourier series, cubic splines)
* the user could even use raw difference in the data as an estimate of the dynamics

$$ \dot{Y} = f(Y,\Omega) $$

### Parameters of process model

```R
K_p   = 10                                            # number of models to fit
W_p   = rep(10,N)                                     # number of neurons in single layer perceptron (SLP)
N_p   = 2 * W_p * (2+N)                               # number of parameters in process model
sd1_p = 0.1                                           # standard deviation of model likelihood
sd2_p = list(c(rep(1.0,N_p[1]/2),rep(.15,N_p[1]/2)),  # standard deviation of prior distributions (second half concerns nonlinear functions)
             c(rep(1.0,N_p[2]/2),rep(.01,N_p[2]/2)),
             c(rep(1.0,N_p[3]/2),rep(.075,N_p[3]/2)))
```

### Train process model

```R
model_p    = trainModel_p(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p)
Yhat_p     = model_p$Yhat_p     
ddx.Yhat_p = model_p$ddx.Yhat_p 
Geber_p    = model_p$Geber_p   
Omega_p    = model_p$Omega_p   
```

### Visualise process model

```R
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
dev.off()
```

### Store results 

```R
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
```

