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
##            - created v0_2
##            - simplified model training

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
for(i in 1:N)
{
	## iterator
	message(paste("fitting: ",i,"/",N,sep=""))

	## best regularisation
	# k = which.max(crossVal_p[[i]][,"logLik_t"])

	## fit model
	Omega_p_i = NULL
	for(k in 1:K_p)
	{   
		## dataloader
		source("m5_loadData_p.r")

		## fit
	    Omega_0   = rnorm(N_p,0,sd2_p[i])
	    Omega_f   = argmax.logPost_p(X_,Y_[,i],Omega_0,sd1_p,sd2_p[i])
	    logPost_0 =        logPost_p(X_,Y_[,i],Omega_0,sd1_p,sd2_p[i])
	    logPost_f =        logPost_p(X_,Y_[,i],Omega_f,sd1_p,sd2_p[i])
	    Omega_p_i  = rbind(Omega_p_i,Omega_f)

		## update
	    message(paste(k,"/",K_p,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
	}
	Omega_p[[i]] = Omega_p_i
			
	## predictions
	Yhat_p[[i]]         =     predict.model_p(X_,Omega_p[[i]])
	ddx.Yhat_p[[i]]     = ddx.predict.model_p(X_,Omega_p[[i]])
	Geber_p[[i]]        = lapply(ddx.Yhat_p[[i]],function(x) x * ddt.X)
	names(Geber_p[[i]]) = c("E.Geber_p","q05.Geber_p","q95.Geber_p")
			
	## de-standardise predictions and store
	Yhat_p[[i]][["E.Yhat_p"]]           = mean_y[i] + sd_y[i] * Yhat_p[[i]][["E.Yhat_p"]]
	Yhat_p[[i]][["q05.Yhat_p"]]         = mean_y[i] + sd_y[i] * Yhat_p[[i]][["q05.Yhat_p"]]
	Yhat_p[[i]][["q95.Yhat_p"]]         = mean_y[i] + sd_y[i] * Yhat_p[[i]][["q95.Yhat_p"]]
	ddx.Yhat_p[[i]][["E.ddx.Yhat_p"]]   = sd_y[i] * t(1/sd_x * t(ddx.Yhat_p[[i]][["E.ddx.Yhat_p"]]))
	ddx.Yhat_p[[i]][["q05.ddx.Yhat_p"]] = sd_y[i] * t(1/sd_x * t(ddx.Yhat_p[[i]][["q05.ddx.Yhat_p"]]))
	ddx.Yhat_p[[i]][["q95.ddx.Yhat_p"]] = sd_y[i] * t(1/sd_x * t(ddx.Yhat_p[[i]][["q95.ddx.Yhat_p"]]))
	Geber_p[[i]][["E.Geber_p"]]         = sd_y[i] * t(1/sd_x * t(Geber_p[[i]][["E.Geber_p"]]))
	Geber_p[[i]][["q05.Geber_p"]]       = sd_y[i] * t(1/sd_x * t(Geber_p[[i]][["q05.Geber_p"]]))
	Geber_p[[i]][["q95.Geber_p"]]       = sd_y[i] * t(1/sd_x * t(Geber_p[[i]][["q95.Geber_p"]]))
}

## store results
save(Yhat_p    ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p   ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p   ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
save(crossVal_p,file=paste(pathToOut,"/","crossVal_p.RData",sep=""))

#
###
