######################
## m0_mainRStudio.r ##
######################

## goal: 
## - analyse NODE to derive interactions and contributions between variables (sensitivities, Geber method)
## - use gradient matching (GM) to fit neural ordinary differential equation model (NODE) to time series 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

## notes:
## - any functions in the supporting function repository "f_NODE_GM_RStudio.r" can be customise in this script
## - for instance the user could specify a custom process model (instead of the single layer perceptron that we use as default)

#
###

##############
## INITIATE ##
##############

## goal: load data, functions

## load data
TS = read.table("data/TS_1.csv",sep=";",header=T)
for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.005)] = 0.005}

## load NODE functions
source("f_NODE_GM_Rstudio.r")

## make output directory
pathToOut = "out/"
system(paste("mkdir",pathToOut))

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model (i.e. interpolate each variable in time series and compute interpolated temporal derivative to approximate temporal dynamics)

## parameters of observation model
K_o     = 100                # number of ensemble elements 
W_o     = c(100,100,100)     # number of neurons in observation model, by default a single layer perceptron (equivalent to number of elements in Fourier series)
N_o     = W_o*3              # total number of parameters in observation model
rho     = 1                  # proportion of best samples to reject (to refine quality of fit if necessary)
alpha_i = 1                  # upsampling interpolation factor (2 double the number of points in interpolated time series)

## train observation model
model_o     = trainModel_o(TS,alpha_i,N_o,K_o,rho)
Yhat_o      = model_o$Yhat_o
ddt.Yhat_o  = model_o$ddt.Yhat_o

## visualise observation model fit
pdf(paste(pathToOut,"/fig_predictions_o.pdf",sep=""))
plotModel_o(TS,alpha_i,Yhat_o,ddt.Yhat_o)
dev.off()

## save results
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model (i.e. explain the per-capita growth rate of the populations calculated as 1/Y*dY/dt as a function of the states Y(t))

## notes: 
## - the user could use state interpolations and interpolated dynamics obtained via other methods (e.g. Fourier series, cubic splines)
## - the user could even use raw difference in the data as an estimate of the dynamics

## parameters of process model
K_p   = 10                                            # number of models to fit
W_p   = rep(10,N)                                     # number of neurons in single layer perceptron (SLP)
N_p   = 2 * W_p * (2+N)                               # number of parameters in process model
sd1_p = 0.1                                           # standard deviation of model likelihood
sd2_p = list(c(rep(1.0,N_p[1]/2),rep(.15,N_p[1]/2)),  # standard deviation of prior distributions (second half concerns nonlinear functions)
             c(rep(1.0,N_p[2]/2),rep(.01,N_p[2]/2)),
             c(rep(1.0,N_p[3]/2),rep(.075,N_p[3]/2)))

## train process model
model_p    = trainModel_p(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p)
Yhat_p     = model_p$Yhat_p     
ddx.Yhat_p = model_p$ddx.Yhat_p 
Geber_p    = model_p$Geber_p   
Omega_p    = model_p$Omega_p   

## visualise process model
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
dev.off()

## store results 
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))

#
###

####################################
## CROSS-VALIDATION PROCESS MODEL ##
####################################

## goal: perform cross-validation of process model by predicting the second half of the time series (of the interpolated p.c. growth rates)
##       the cross-validation is performed to identify the best value for regularisation on the nonlinear part of the model

## parameters for cross-validation
K_p             = 3                                      # number of models to fit per folds and regularisation parameter
folds           = list(c(1/2,1))                         # proportion of the data that should be considered for training (e.g. c(1/2,1) means that the second half is for testing) 
# folds           = list(c(0,1/3),c(1/3,2/3),c(2/3,3/3)) # alternative folds
crossValParVect = seq(0.01,0.5,0.01)

## run cross-validation
resultsCrossVal_p = crossVal_p(TS,alpha_i,Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,folds,crossValParVect)

## visualise
pdf(paste(pathToOut,"/fig_crossVal_p.pdf",sep=""))
plotCrossVal_p(resultsCrossVal_p)
dev.off()

## store results
save(resultsCrossVal_p,file=paste(pathToOut,"/","crossVal_p.RData"   ,sep=""))

#
###
