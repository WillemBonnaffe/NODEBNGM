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

	## best regularisation
	# k = which.max(crossVal_p[[i]][,"logLik_t"])

	## fit model
	Omega_p_i = NULL
	for(k in 1:K_p)
	{   

		## dataloader 
		source("m5_loadData_p.r")

		## fit
	    Omega_0      = rnorm(N_p,0,sd2_p[i])
        Yhat         = function(X,Omega) f_p.eval(X,Omega,w[i])
        ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega,w[i])
	    # Omega_f   = argmax.logPost(X_,Y_[,i],Yhat,ddOmega.f_p.eval,Omega_0,sd1_p,sd2_p[i])
	    Omega_f      = argmax.logMarPost(X_,Y_[,i],Yhat,ddOmega.Yhat,Omega_0)
	    Omega_p_i    = rbind(Omega_p_i,Omega_f)

		## update
	    # logPost_0  = logPost(X_,Y_[,i],Yhat,Omega_0,sd1_p,sd2_p[i])
	    # logPost_f  = logPost(X_,Y_[,i],Yhat,Omega_f,sd1_p,sd2_p[i])
        logPost_0    = logMarPost(X_,Y_[,i],Yhat,Omega_0)
	    logPost_f    = logMarPost(X_,Y_[,i],Yhat,Omega_f)
	    message(paste(k,"/",K_p,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
	}

	## store
	Omega_p[[i]] = Omega_p_i
			
	## predictions
    Yhat_p[[i]]     = t(apply(Omega_p[[i]],1,function(x) mean_y[i] + sd_y[i] * f_p.eval(X_,x,w[i])))
	ddx.Yhat_p[[i]] = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * ddx.f_p.eval(X_,x,w[i])))
	Geber_p[[i]]    = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * t(E.ddt.X_p) * ddx.f_p.eval(X_,x,w[i])))
}

## store results
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))

#
###
