#######################
## m6_trainModel_p.r ##
#######################

## goal: apply process model to analyse interaction between variables driving temporal dynamics of system 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0
## 14-04-2022 - created v0_1
##            - loading different observation data at every iteration
## 15-04-2022 - simplified code
##            - evaluate predictions for all Omega_p_ik

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("m5_loadModel_p.r")

#
###

#########################
## TRAIN PROCESS MODEL ##
#########################

## goal: fit process model

## for each variable 
Omega_p    = list()
Yhat_p     = list()
ddx.Yhat_p = list()
Geber_p    = list()
crossVal_p = list()
for(i in 1:N)
{
	## iterator
	message(paste("fitting: ",i,"/",N,sep=""))
			 	
	## cross validation
	crossVal_i      = NULL
	Omega_p[[i]]    = list()
	Yhat_p[[i]]     = list()
	ddx.Yhat_p[[i]] = list()
	Geber_p[[i]]    = list()
	for(k in 1:length(sd2_p)) # for each regularisation param
	{
		## iterator
		message(paste("crossval: ",k,"/",length(sd2_p),sep=""))

		## multiple folds
		crossVal_ik          = NULL
		Omega_p[[i]][[k]]    = list()
		Yhat_p[[i]][[k]]     = list()
		ddx.Yhat_p[[i]][[k]] = list()
		Geber_p[[i]][[k]]    = list()
		for(u in 1:length(folds)) # for each fold in the data
		{
			## iterator
			message(paste("fold: ",u,"/",length(folds),sep=""))

			## dataloader
			source("m5_loadData_p.r")
			X_lu   = X_l[[u]]
			Y_lu   = Y_l[[u]][,i]
			X_tu   = X_t[[u]]
			Y_tu   = Y_t[[u]][,i]
			sd_y   = sd_y[i]
			mean_y = mean_y[i]

			## fit model_p
			Omega_p[[i]][[k]][[u]] = fit.model_p(X_lu,Y_lu,N_p,sd1_p,sd2_p[[k]],K_p)
			
			## model performance
			logLik_l = apply(Omega_p[[i]][[k]][[u]],1,function(x)logLik_p(X_lu,Y_lu,x,sd1_p))
			logLik_t = apply(Omega_p[[i]][[k]][[u]],1,function(x)logLik_p(X_tu,Y_tu,x,sd1_p))
					
			## cross validation matrix
			crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))

			## predictions
			Yhat_p[[i]][[k]][[u]]         =     predict.model_p(X_,Omega_p[[i]][[k]][[u]])
			ddx.Yhat_p[[i]][[k]][[u]]     = ddx.predict.model_p(X_,Omega_p[[i]][[k]][[u]])
			Geber_p[[i]][[k]][[u]]        = lapply(ddx.Yhat_p[[i]],function(x) x * ddt.X)
			names(Geber_p[[i]][[k]][[u]]) = c("E.Geber_p","q05.Geber_p","q95.Geber_p")
			
			## de-standardise predictions and store
			Yhat_p[[i]][[k]][[u]][["E.Yhat_p"]]           = mean_y + sd_y * Yhat_p[[i]][[k]][[u]][["E.Yhat_p"]]
			Yhat_p[[i]][[k]][[u]][["q05.Yhat_p"]]         = mean_y + sd_y * Yhat_p[[i]][[k]][[u]][["q05.Yhat_p"]]
			Yhat_p[[i]][[k]][[u]][["q95.Yhat_p"]]         = mean_y + sd_y * Yhat_p[[i]][[k]][[u]][["q95.Yhat_p"]]
			ddx.Yhat_p[[i]][[k]][[u]][["E.ddx.Yhat_p"]]   = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][[k]][[u]][["E.ddx.Yhat_p"]]))
			ddx.Yhat_p[[i]][[k]][[u]][["q05.ddx.Yhat_p"]] = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][[k]][[u]][["q05.ddx.Yhat_p"]]))
			ddx.Yhat_p[[i]][[k]][[u]][["q95.ddx.Yhat_p"]] = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][[k]][[u]][["q95.ddx.Yhat_p"]]))
			Geber_p[[i]][[k]][[u]][["E.Geber_p"]]         = sd_y * t(1/sd_x * t(Geber_p[[i]][[k]][[u]][["E.Geber_p"]]))
			Geber_p[[i]][[k]][[u]][["q05.Geber_p"]]       = sd_y * t(1/sd_x * t(Geber_p[[i]][[k]][[u]][["q05.Geber_p"]]))
			Geber_p[[i]][[k]][[u]][["q95.Geber_p"]]       = sd_y * t(1/sd_x * t(Geber_p[[i]][[k]][[u]][["q95.Geber_p"]]))

		}

		## store 
		E.crossVal_ik  = apply(crossVal_ik,2,mean)
		sd.crossVal_ik = apply(crossVal_ik,2,sd)
		crossVal_i     = rbind(crossVal_i,c(sd2_p[[k]],E.crossVal_ik,sd.crossVal_ik))
	}

	## store
	crossVal_p[[i]] = crossVal_i
	colnames(crossVal_p[[i]]) = c("sd","logLik_l","logLik_t","sd.logLik_l","sd.logLik_t")

	## best regularisation
	k        = which.max(crossVal_p[[i]][,"logLik_t"])
	Omega_pi = NULL
	for(u in 1:length(folds))
	{
		 Omega_pi = rbind(Omega_pi,Omega_p[[i]][[k]][[u]])
	}	
}

## store results
save(Yhat_p    ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p   ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p   ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
save(crossVal_p,file=paste(pathToOut,"/","crossVal_p.RData",sep=""))

#
###
