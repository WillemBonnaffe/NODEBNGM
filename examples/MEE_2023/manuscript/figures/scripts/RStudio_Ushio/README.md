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

### Parameters of the observation process model

The second input is the parameters of the observation model (i.e. the neural networks that interpolates the time series) and process model (the NODEs that approximate the per-capita growth rate based on the interpolated variables).
See below for instructions on how to specify the parameters and run the code.


## Outputs

### Interpolated time series and dynamics

### Effects and contributions of variables to the dynamics of other variables

### Dynamical interaction networks


## Installation

The approach can be used simply by installing R (v4.0.2 or later) and loading the NODEBNGM function library

``` R
source("f_NODE_GM_Rstudio.r")
```


## Preparing the data

### Loading data

``` R
TS = read.table("data/TS.csv",sep=",",header=T)
```

### Select variables of interest

``` R
selected_time_steps = 50:150
selected_columns  = c(
  "time_step",
  # "surf.t",
  "bot.t",
  "Aurelia.sp",
  # "Engraulis.japonicus", #
  # "Plotosus.lineatus", #
  "Sebastes.inermis",
  "Trachurus.japonicus",
  "Girella.punctata",
  "Pseudolabrus.sieboldi",
  "Halichoeres.poecilopterus",
  "Halichoeres.tenuispinnis",
  # "Chaenogobius.gulosus", #
  "Pterogobius.zonoleucus",
  "Tridentiger.trigonocephalus",
  # "Siganus.fuscescens", #
  "Sphyraena.pinguis", #
  "Rudarius.ercodes"           
)
TS = TS[selected_time_steps,]
TS = TS[,selected_columns]
```

### Format column names

```R
column_names =  c("time_step",
                    "bot.t",
                    "Aurel.sp",
                    "S.inerm.",
                    "T.japon.",
                    "G.punct.",
                    "P.siebo.",
                    "H.poeci.",
                    "H.tenui.",
                    "P.zonol.",
                    "T.trigo.",
                    "S.pingu.",
                    "R.ercod."
                  )
colnames(TS) = column_names
```

### Normalising time series

``` R
TS[,-1] = apply(TS[,-1],2,function(x)(x-min(x))/(max(x)-min(x))*10) # normalise time series
for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.005)] = 0.005} # set 0s to small value to avoid NAs
```

### Visualise the formatted time series

``` R
par(mfrow=c(3,4))
for(i in 2:ncol(TS))
{
  plot(TS[,1],TS[,i],type="l",xlab="Time step",ylab="Relative density",bty="n",main=colnames(TS)[i])
}
par(mfrow=c(1,1))
```

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_time_series.png)


## Fitting observation model

Fitting the observation model corresponds in interpolating the variables in the time series to get an approximation of the states and dynamics of each variable.

### Parameters of the observation model

``` R
N       = ncol(TS) - 1
K_o     = 100                # number of ensemble elements
W_o     = rep(30,N)          # number of neurons in observation model, by default a single layer perceptron (equivalent to number of elements in Fourier series)
N_o     = W_o*3              # total number of parameters in observation model
rho     = 1                  # proportion of best samples to reject (to refine quality of fit if necessary)
alpha_i = 1                  # upsampling interpolation factor (2 double the number of points in interpolated time series)
```

### Training of the observation model

``` R
model_o     = trainModel_o(TS,alpha_i,N_o,K_o,rho)
Yhat_o      = model_o$Yhat_o
ddt.Yhat_o  = model_o$ddt.Yhat_o
Omega_o     = model_o$Omega_o
```

### Visualising the fit

``` R
plotModel_o(TS,alpha_i,Yhat_o,ddt.Yhat_o)
```

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_predictions_o.png)

### Storing the results

``` R
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))
```


## Fitting process model

Fit process model (i.e. explain the per-capita growth rate of the populations calculated as $1/Y*dY/dt$ as a function of the states Y(t))

### Load observation model results (optional)

``` R
load(file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
load(file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
load(file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))
```

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
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
```

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_predictions_p.png)


### Store results 

```R
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
```


## Analysing results

This section describes how to run the code to obtain effects and contributions from the results obtained from the observation and process model.

### Load process model results

``` R
load(paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
load(paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
load(paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
```

### Compute Jacobian matrix

```R
## compute Jacobian and contribution matrix
MSq = function(x) mean(x^2)
prop = function(x) x/sum(x)
J = t(matrix(unlist(lapply(ddx.Yhat_p,function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,mean))),ncol=ncol(TS)-1)) ## average across samples then average across time steps
C = t(matrix(unlist(lapply(Geber_p,   function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,MSq))),ncol=ncol(TS)-1)) ## average across samples then take mean square across time steps
C = t(apply(C,1,prop))
```

### Formatting the Jacobian matrix

``` R
## remove effects on bot
J[1,] = 0
C[1,] = 0

## thresholding
J = J*(C>0.1)
C = C*(C>0.1)
```

### Visualise the dynamical interaction plot

``` R
##
## DYNAMICAL INTERACTION PLOT (v1)
.plot.DIN(J,C,colnames(TS)[-1])
```

![alt text](https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/manuscript/figures/scripts/RStudio_Ushio/out/fig_DIN_v1.png)

## Notes 
* the user could use state interpolations and interpolated dynamics obtained via other methods (e.g. Fourier series, cubic splines)
* the user could even use raw difference in the data as an estimate of the dynamics
