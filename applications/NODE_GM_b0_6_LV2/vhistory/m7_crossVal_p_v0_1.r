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

## parameters
K_p             = 30
folds           = list(c(0,1/3),c(1/3,2/3),c(2/3,3/3))
crossValParVect = seq(0.1,0.5,0.1)

#
###

#########################
## TRAIN PROCESS MODEL ##
#########################

## goal: crossval process model

## for each variable 
Omega_p    = list()
crossVal_p = list()
for(i in 1:N)
{
	## iterator
	message(paste("fitting: ",i,"/",N,sep=""))
			 	
	## cross validation
	crossVal_i   = NULL
	Omega_p[[i]] = list()
	for(k in 1:length(crossValParVect)) # for each regularisation param
	{
		## iterator
		message(paste("crossval: ",k,"/",length(crossValParVect),sep=""))

		## multiple folds
		crossVal_ik          = NULL
		Omega_p[[i]][[k]]    = list()
		for(u in 1:length(folds)) # for each fold in the data
		{
			## iterator
			message(paste("fold: ",u,"/",length(folds),sep=""))

			## fit model
			Omega_p_iku = NULL
			for(m in 1:K_p)
			{   
				## dataloader
				source("m5_loadData_p.r")

				## split training/test
			    s   = round(folds[[u]][1]*n+1):round(folds[[u]][2]*n)
				X_l = X_[-s,]
				Y_l = Y_[-s,]
				X_t = X_[s,]
				Y_t = Y_[s,]

			    ## fit
			    Omega_0     = rnorm(N_p,0,crossValParVect[k])
			    Omega_f     = argmax.logPost_p(X_l,Y_l[,i],Omega_0,sd1_p,crossValParVect[k])
			    Omega_p_iku = rbind(Omega_p_iku,Omega_f)

				## update
				logPost_0 = logPost_p(X_l,Y_l[,i],Omega_0,sd1_p,crossValParVect[k])
				logPost_f = logPost_p(X_l,Y_l[,i],Omega_f,sd1_p,crossValParVect[k])
				# message(paste(m,"/",K_p,"\t",
        		# 		format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
        		# 		format(round(logPost_f,2),nsmall=2),sep=""))

				## model performance
				logLik_l = logLik_p(X_l,Y_l[,i],Omega_f,sd1_p)
				logLik_t = logLik_p(X_t,Y_t[,i],Omega_f,sd1_p)						
				# message(paste(m,"/",K_p,"\t",
        		# 		format(round(logLik_l,2),nsmall=2),"\t","","\t",
        		# 		format(round(logLik_t,2),nsmall=2),sep=""))

				## cross validation matrix
				crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))

			}       
			Omega_p[[i]][[k]][[u]] = Omega_p_iku
		}

		## store 
		E.crossVal_ik  = apply(crossVal_ik,2,mean)
		sd.crossVal_ik = apply(crossVal_ik,2,sd)
		crossVal_i     = rbind(crossVal_i,c(crossValParVect[k],E.crossVal_ik,sd.crossVal_ik))
		message(paste("logLik l vs t: ",
        		format(round(E.crossVal_ik[1],2),nsmall=2),"\t",
        		format(round(E.crossVal_ik[2],2),nsmall=2),sep=""))

	}

	## store
	crossVal_p[[i]] = crossVal_i
	colnames(crossVal_p[[i]]) = c("sd","logLik_l","logLik_t","sd.logLik_l","sd.logLik_t")
}

## store results
save(Omega_p   ,file=paste(pathToOut,"/","crossVal_p-Omega_p.RData"   ,sep=""))
save(crossVal_p,file=paste(pathToOut,"/","crossVal_p.RData",sep=""))

#
###
