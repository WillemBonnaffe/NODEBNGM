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
## 25-04-2022 - created v0_2
##            - not relying on function fit anymore
##            - simplified code by transforming ensemble rather than each derived quantity
## 03-06-2022 - created v0_0 from m3_trainModel.r
##            - to implement AMC
##            - implemented DEMC

##############
## INITIATE ##
##############

## goal: initiate the observation model

##
library("Rcpp")
sourceCpp("cpp/DEMCpp_v0.1.cpp")

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

	## get ensemble
	Omega_o_i = NULL
	for(k in 1:1)
	{   
		## fit
	    Omega_0 = rnorm(N_o,0,sd2_o[i])
        dTarget = function(x) logMarPost(t_,Y_[,i],f_o.eval,x,1/W_o)
        chain   = DEMCpp(list("dTarget" = dTarget,
        					   "Theta_0" = Omega_0, 
        					   "epsilon" = 0.001, 
        					   "nIt"     = 100000))$chainList
        Omega_o_i = chain[,-1]
	}

	## store ensemble
	Omega_o[[i]] = Omega_o_i

	## compute predictions
	Yhat_o[[i]]     = t(apply(Omega_o[[i]],1,function(x)     f_o.eval(nt_,x))) 
	ddt.Yhat_o[[i]] = t(apply(Omega_o[[i]],1,function(x) ddt.f_o.eval(nt_,x)))

	## de-standardise
	Yhat_o[[i]]     = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]])
	ddt.Yhat_o[[i]] = 1/dt * sd_y[i] * Yhat_o[[i]] * ddt.Yhat_o[[i]]
 }

## store results
names(Yhat_o)     = colnames(TS[,-1])
names(ddt.Yhat_o) = colnames(TS[,-1])
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###
