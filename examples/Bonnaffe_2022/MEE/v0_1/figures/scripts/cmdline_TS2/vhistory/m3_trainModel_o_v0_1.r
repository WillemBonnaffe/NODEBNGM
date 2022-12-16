#######################
## m3_trainModel_o.r ##
#######################

## goal: fit observation model to interpolate and estimate temporal derivative of time series data 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0
## 14-04-2022 - created v0_1
##            - added transformation of ensemble
## 15-04-2022 - added dataloader
##            - simplified code

##############
## INITIATE ##
##############

## goal: initiate the observation model

## load data and model 
source("m1_loadData_o.r")
source("m2_loadModel_o.r")

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model 

## interpolate each column in the time series
Omega_o    = list()
Yhat_o     = list()
ddt.Yhat_o = list()
for (i in 1:N)
{
 	## iterator
 	message(paste("fitting: ",i,"/",ncol(Y_o),sep=""))

	## fit model_o (get ensemble of argmax)
	Omega_o[[i]] = fit.model_o(t_,Y_[,i],N_o,sd1_o,sd2_o,K_o,logMar=T)

	## compute predictions
	Yhat_o[[i]]     =     predict.model_o(nt_,Omega_o[[i]])
	ddt.Yhat_o[[i]] = ddt.predict.model_o(nt_,Omega_o[[i]])

	## de-standardise (for logged variables) (remove exp in Yhat_o and Yhat_o in ddt.Yhat_o for non-logged variables)
	Yhat_o[[i]][["MaP.Yhat_o"]]         = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]][["MaP.Yhat_o"]])
	Yhat_o[[i]][["E.Yhat_o"]]           = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]][["E.Yhat_o"]])
	Yhat_o[[i]][["q05.Yhat_o"]]         = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]][["q05.Yhat_o"]])
	Yhat_o[[i]][["q95.Yhat_o"]]         = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]][["q95.Yhat_o"]])
	Yhat_o[[i]][["Yhat_o"]]             = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]][["Yhat_o"]])
	ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]] = 1/dt * sd_y[i] * Yhat_o[[i]][["MaP.Yhat_o"]] * ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]]
	ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]   = 1/dt * sd_y[i] * Yhat_o[[i]][["E.Yhat_o"]]   * ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]
	ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]] = 1/dt * sd_y[i] * Yhat_o[[i]][["E.Yhat_o"]]   * ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]]
	ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]] = 1/dt * sd_y[i] * Yhat_o[[i]][["E.Yhat_o"]]   * ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]]
	ddt.Yhat_o[[i]][["ddt.Yhat_o"]]     = 1/dt * sd_y[i] * Yhat_o[[i]][["Yhat_o"]]     * ddt.Yhat_o[[i]][["ddt.Yhat_o"]]
 }

## store results
names(Yhat_o)     = colnames(TS[,-1])
names(ddt.Yhat_o) = colnames(TS[,-1])
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###
