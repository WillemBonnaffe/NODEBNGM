#######################
## m7_trainModel_p.r ##
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
## 25-04-2022 - created v0_3
##            - re-introduced marginal posterior
## 30-05-2022 - created v0_4
##            - implemented standard regularisation again with specific priors for each variables
##            - this is to control to tease apart the importance of nonlinearity in individual variables
##            - implemented marginal posterior again as greatly improves fitting
##            - implemented variable specific network size (cross validation performed on network size instead of "slider")
##            - removed slider given that hyperparameter is now size of network

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("m6_loadModel_p.r")

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

	## fit model
	Omega_p_i = NULL
	for(k in 1:K_p)
	{   

		## dataloader 
		source("m5_loadData_p.r")

		## fit
        N_p          = 2 * W_p[i] * (2+N)
	    Omega_0      = rnorm(N_p,0,0.1)
        Yhat         = function(X,Omega) f_p.eval(X,Omega)
        ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega)
	    Omega_f      = argmax.logMarPost(X_,Y_[,i],Yhat,ddOmega.Yhat,Omega_0,1/W_p[i])
	    Omega_p_i    = rbind(Omega_p_i,Omega_f)

		## update
	    logPost_0    = logMarLik(X_,Y_[,i],Yhat,Omega_0)
	    logPost_f    = logMarLik(X_,Y_[,i],Yhat,Omega_f)
	    message(paste(k,"/",K_p,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
	}

	## store
	Omega_p[[i]] = Omega_p_i
			
	## predictions
    Yhat_p[[i]]     = t(apply(Omega_p[[i]],1,function(x) mean_y[i] + sd_y[i] * f_p.eval(X_,x)))
	ddx.Yhat_p[[i]] = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * ddx.f_p.eval(X_,x)))
	Geber_p[[i]]    = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * t(E.ddt.X_p) * ddx.f_p.eval(X_,x)))
}

## store results
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))

#
###
