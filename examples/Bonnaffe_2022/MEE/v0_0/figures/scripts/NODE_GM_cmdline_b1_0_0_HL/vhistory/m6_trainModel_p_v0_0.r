#######################
## m6_trainModel_p.r ##
#######################

## goal: apply process model to analyse interaction between variables driving temporal dynamics of system 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0

#########################
## TRAIN PROCESS MODEL ##
#########################

## goal: fit process model

## iterator 
Omega_p    = list()
Yhat_p     = list()
ddx.Yhat_p = list()
Geber_p    = list()
crossVal_p = list()
for(i in 1:ncol(Y_p))
{
	## iterator
	message(paste("fitting: ",i,"/",ncol(Y_p),sep=""))

	## variables
	X     = X_p 
	ddt.X = ddt.X_p
	Y     = Y_p[,i]
	
	## standardise predictive variables
	X_     = X
	mean_x = apply(X_,2,mean)
	sd_x   = apply(X_,2,sd)
	X_     = t((t(X_)-mean_x)/sd_x)
	
	## standardise response variable
	Y_     = Y
	mean_y = mean(Y_)
	sd_y   = sd(Y_)
	Y_     = (Y_-mean_y)/sd_y
			 	
	## cross validation
	crossVal_i   = NULL
	Omega_p[[i]] = list()
	for(k in 1:length(sd2_p)) # for each regularisation param
	{
		## iterator
		message(paste("crossval: ",k,"/",length(sd2_p),sep=""))

		## multiple folds
		crossVal_ik       = NULL
		Omega_p[[i]][[k]] = list()
		for(u in 1:length(folds)) # for each fold in the data
		{
			## iterator
			message(paste("fold: ",u,"/",length(folds),sep=""))

			## split training/test
			s   = round(folds[[u]][1]*length(Y_)+1):round(folds[[u]][2]*length(Y_))
			X_l = X_[-s,]
			Y_l = Y_[-s]
			X_t = X_[s,]
			Y_t = Y_[s]

			## fit model_p
			Omega_p[[i]][[k]][[u]] = fit.model_p(X_l,Y_l,N_p,sd1_p,sd2_p[[k]],K_p)
			
			## model performance
			logLik_l = apply(Omega_p[[i]][[k]][[u]],1,function(x)logLik_p(X_l,Y_l,x,sd1_p))
			logLik_t = apply(Omega_p[[i]][[k]][[u]],1,function(x)logLik_p(X_t,Y_t,x,sd1_p))
					
			## cross validation matrix
			crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))
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
	
	## predictions
	Yhat_p[[i]]         = predict.model_p(X_,Omega_pi)
	ddx.Yhat_p[[i]]     = ddx.predict.model_p(X_,Omega_pi)
	Geber_p[[i]]        = lapply(ddx.Yhat_p[[i]],function(x) x * ddt.X)
	names(Geber_p[[i]]) = c("E.Geber_p","q05.Geber_p","q95.Geber_p")
	
	## de-standardise predictions and store
	Yhat_p[[i]][["E.Yhat_p"]]           = mean_y + sd_y * Yhat_p[[i]][["E.Yhat_p"]]
	Yhat_p[[i]][["q05.Yhat_p"]]         = mean_y + sd_y * Yhat_p[[i]][["q05.Yhat_p"]]
	Yhat_p[[i]][["q95.Yhat_p"]]         = mean_y + sd_y * Yhat_p[[i]][["q95.Yhat_p"]]
	ddx.Yhat_p[[i]][["E.ddx.Yhat_p"]]   = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][["E.ddx.Yhat_p"]]))
	ddx.Yhat_p[[i]][["q05.ddx.Yhat_p"]] = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][["q05.ddx.Yhat_p"]]))
	ddx.Yhat_p[[i]][["q95.ddx.Yhat_p"]] = sd_y * t(1/sd_x * t(ddx.Yhat_p[[i]][["q95.ddx.Yhat_p"]]))
	Geber_p[[i]][["E.Geber_p"]]         = sd_y * t(1/sd_x * t(Geber_p[[i]][["E.Geber_p"]]))
	Geber_p[[i]][["q05.Geber_p"]]       = sd_y * t(1/sd_x * t(Geber_p[[i]][["q05.Geber_p"]]))
	Geber_p[[i]][["q95.Geber_p"]]       = sd_y * t(1/sd_x * t(Geber_p[[i]][["q95.Geber_p"]]))
}

## store results
save(Yhat_p    ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p   ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p   ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
save(crossVal_p,file=paste(pathToOut,"/","crossVal_p.RData",sep=""))

#
###
