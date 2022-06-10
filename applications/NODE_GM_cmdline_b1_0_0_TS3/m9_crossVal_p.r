#####################
## m9_crossVal_p.r ##
#####################

## goal: perform cross validation 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-04-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("m6_loadModel_p.r")

## parameters
K_p             = 10
folds           = list(c(1/2,1)) # list(c(0,1/3),c(1/3,2/3),c(2/3,3/3))
crossValParVect = seq(0.01,.3,0.025)

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

                ## regularisation on nonlinear part of model
                sd2_p[[i]][(N_p[i]/2+1):N_p[i]] = rep(crossValParVect[k],N_p[i]/2) 

			    ## fit
			    Omega_0      = rnorm(N_p[i],0,0.01)
                Yhat         = function(X,Omega) f_p.eval(X,Omega)
                ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega)
			    # Omega_f      = argmax.logMarPost(X_l,Y_l[,i],Yhat,ddOmega.Yhat,Omega_0,1/W_p[i])
			    Omega_f      = argmax.logPost(X_l,Y_l[,i],Yhat,ddOmega.Yhat,Omega_0,sd1_p,sd2_p[[i]])
			    Omega_p_iku  = rbind(Omega_p_iku,Omega_f)

				## update
				logPost_0 = logPost(X_l,Y_l[,i],Yhat,Omega_0,sd1_p,sd2_p[[i]])
				logPost_f = logPost(X_l,Y_l[,i],Yhat,Omega_f,sd1_p,sd2_p[[i]])
                # logMarPost_0 = logMarPost(X_l,Y_l[,i],Yhat,Omega_0,1/W_p[i])
                # logMarPost_f = logMarPost(X_l,Y_l[,i],Yhat,Omega_f,1/W_p[i])				
                # message(paste(m,"/",K_p,"\t",
        		# 		format(round(logMarPost_0,2),nsmall=2),"\t","-->","\t",
        		# 		format(round(logMarPost_f,2),nsmall=2),sep=""))

				## model performance
				logLik_l = logLik(X_l,Y_l[,i],Yhat,Omega_f,sd1_p)
				logLik_t = logLik(X_t,Y_t[,i],Yhat,Omega_f,sd1_p)
				# logMarLik_l = logMarLik(X_l,Y_l[,i],Yhat,Omega_f)
				# logMarLik_t = logMarLik(X_t,Y_t[,i],Yhat,Omega_f)
				# message(paste(m,"/",K_p,"\t",
        		# 		format(round(logMarLik_l,2),nsmall=2),"\t","","\t",
        		# 		format(round(logMarLik_t,2),nsmall=2),sep=""))

				## cross validation matrix
				crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))
				# crossVal_ik  = rbind(crossVal_ik,cbind(logMarLik_l,logMarLik_t))

			}       
			Omega_p[[i]][[k]][[u]] = Omega_p_iku
		}

		## store 
		E.crossVal_ik  = apply(crossVal_ik,2,mean)
		sd.crossVal_ik = apply(crossVal_ik,2,sd)
		crossVal_i     = rbind(crossVal_i,c(crossValParVect[k],E.crossVal_ik,sd.crossVal_ik))
		message(paste("logLik l vs t: ",
        		format(round(E.crossVal_ik[1],2),nsmall=2),"\t",
        		format(round(E.crossVal_ik[2],2),nsmall=2),"\n",sep=""))

	}

	## store
	crossVal_p[[i]] = crossVal_i
	colnames(crossVal_p[[i]]) = c("sd","logLik_l","logLik_t","sd.logLik_l","sd.logLik_t")
	# colnames(crossVal_p[[i]]) = c("w","logMarLik_l","logMarLik_t","sd.logMarLik_l","sd.logMarLik_t")
}

## store results
save(Omega_p   ,file=paste(pathToOut,"/","crossVal_p-Omega_p.RData"   ,sep=""))
save(crossVal_p,file=paste(pathToOut,"/","crossVal_p.RData",sep=""))

#
###
