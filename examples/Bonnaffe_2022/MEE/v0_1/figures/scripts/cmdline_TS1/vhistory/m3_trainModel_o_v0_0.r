#######################
## m3_trainModel_o.r ##
#######################

## goal: fit observation model to interpolate and estimate temporal derivative of time series data 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate the observation model

## imports
source("m1_loadData.r")
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

	## predictive and response variable
	t   = X_o 
	Y   = Y_o[,i] 
	nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)

	## standardise time steps
	t_0 = min(t)
	t_f = max(t)
	dt  = diff(t[1:2])
	t_  = 1:length(t)
	nt_ = seq(min(t_),max(t_),(t_[2]-t_[1])/alpha_i)
	
	## standardise data
	Y_     = Y
	if (log[i] == T) Y_ = log(Y_)
	mean_y = mean(Y_)
	sd_y   = sd(Y_)
	Y_     = (Y_-mean_y)/sd_y
	
	## fit model_o (get ensemble of argmax)
	Omega_o[[i]] = fit.model_o(t_,Y_,N_o,sd1_o,sd2_o,K_o,logMar=F)

	## compute predictions
	Yhat_o[[i]]     =     predict.model_o(nt_,Omega_o[[i]])
	ddt.Yhat_o[[i]] = ddt.predict.model_o(nt_,Omega_o[[i]])

	## de-standardise (for logged variables)
	if(log[i] == T)
	{
		Yhat_o[[i]][["MaP.Yhat_o"]]         = exp(mean_y + sd_y * Yhat_o[[i]][["MaP.Yhat_o"]])
		Yhat_o[[i]][["E.Yhat_o"]]           = exp(mean_y + sd_y * Yhat_o[[i]][["E.Yhat_o"]])
		Yhat_o[[i]][["q05.Yhat_o"]]         = exp(mean_y + sd_y * Yhat_o[[i]][["q05.Yhat_o"]])
		Yhat_o[[i]][["q95.Yhat_o"]]         = exp(mean_y + sd_y * Yhat_o[[i]][["q95.Yhat_o"]])
		ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]] = 1/dt * sd_y * Yhat_o[[i]][["MaP.Yhat_o"]] * ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]]
		ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]   = 1/dt * sd_y * Yhat_o[[i]][["E.Yhat_o"]]   * ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]
		ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]] = 1/dt * sd_y * Yhat_o[[i]][["E.Yhat_o"]] * ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]]
		ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]] = 1/dt * sd_y * Yhat_o[[i]][["E.Yhat_o"]] * ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]]
	}
	else
	{	
	 	## de-standardise (for natural variables)
	 	Yhat_o[[i]][["MaP.Yhat_o"]]         = mean_y + sd_y *     Yhat_o[[i]][["MaP.Yhat_o"]]
	 	Yhat_o[[i]][["E.Yhat_o"]]           = mean_y + sd_y *     Yhat_o[[i]][["E.Yhat_o"]]  
	 	Yhat_o[[i]][["q05.Yhat_o"]]         = mean_y + sd_y *     Yhat_o[[i]][["q05.Yhat_o"]]
	 	Yhat_o[[i]][["q95.Yhat_o"]]         = mean_y + sd_y *     Yhat_o[[i]][["q95.Yhat_o"]]
	 	ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]] =   1/dt * sd_y * ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]]
	 	ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]   =   1/dt * sd_y * ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]]
	 	ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]] =   1/dt * sd_y * ddt.Yhat_o[[i]][["q05.ddt.Yhat_o"]]
	 	ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]] =   1/dt * sd_y * ddt.Yhat_o[[i]][["q95.ddt.Yhat_o"]]
	}
}

## store results
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###
