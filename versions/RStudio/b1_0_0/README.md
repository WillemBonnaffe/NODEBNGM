############
## README ##
############

## update log:
## 28-03-2022 - create b0_6

#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model (i.e. explain the per-capita growth rate of the populations calculated as 1/Y*dY/dt as a function of the states Y(t))

## notes: 
## - the user could use state interpolations and interpolated dynamics obtained via other methods (e.g. Fourier series, cubic splines)
## - the user could even use raw difference in the data as an estimate of the dynamics

## parameters of process model
'''
K_p   = 10                                            # number of models to fit
W_p   = rep(10,N)                                     # number of neurons in single layer perceptron (SLP)
N_p   = 2 * W_p * (2+N)                               # number of parameters in process model
sd1_p = 0.1                                           # standard deviation of model likelihood
sd2_p = list(c(rep(1.0,N_p[1]/2),rep(.15,N_p[1]/2)),  # standard deviation of prior distributions (second half concerns nonlinear functions)
             c(rep(1.0,N_p[2]/2),rep(.01,N_p[2]/2)),
             c(rep(1.0,N_p[3]/2),rep(.075,N_p[3]/2)))
'''

## train process model

'''
model_p    = trainModel_p(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p)
Yhat_p     = model_p$Yhat_p     
ddx.Yhat_p = model_p$ddx.Yhat_p 
Geber_p    = model_p$Geber_p   
Omega_p    = model_p$Omega_p   
'''

## visualise process model
'''
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
dev.off()

## store results 
'''
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
'''

#
###
