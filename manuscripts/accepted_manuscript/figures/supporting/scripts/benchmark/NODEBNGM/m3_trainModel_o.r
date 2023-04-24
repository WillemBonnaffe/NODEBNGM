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
## 08-06-2022 - created v0_3
##            - implemented rejection of samples too close to initial conditions

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
logPost_o  = list()
runtimes   = NULL
for (i in 1:N)
{
 	## iterator
 	message(paste("fitting: ",i,"/",ncol(Y_o),sep=""))

	## get ensemble
	Omega_o_i = NULL
    logPost_o_i = NULL
	for(k in 1:K_o)
	{   
        check = F
        while(check == F)
        {
		    ## fit
	        Omega_0   = rnorm(N_o[i],0,sd2_o[i])
	        runtime = system.time({
	            Omega_f   = argmax.logMarPost(t_,Y_[,i],f_o.eval,ddOmega.f_o.eval,Omega_0,1/W_o[i])
	        })[3]

		    ## update
            logPost_0 = logMarPost(t_,Y_[,i],f_o.eval,Omega_0,1/W_o[i])
	        logPost_f = logMarPost(t_,Y_[,i],f_o.eval,Omega_f,1/W_o[i])
	        message(paste(k,"/",K_o,"\t",
	                format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	                format(round(logPost_f,2),nsmall=2),sep=""))
            
            ## reject or keep sample
            check = (logPost_f >= logPost_0 + 10)
            if(check == T)
            {
	   	        Omega_o_i = rbind(Omega_o_i,Omega_f)
                # logPost_o_i = c(logPost_o_i,logPost_f)
                logPost_o_i = c(logPost_o_i,logMarLik(t_,Y_[,i],f_o.eval,Omega_f))
            }
            
            ## store runtimes
            runtimes = c(runtimes, runtime)
        }
	}

    ## sort by best models
    s = rev(order(logPost_o_i))[1:5]
    Omega_o_i = Omega_o_i[s,]

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
write(runtimes, file=paste(pathToOut,"/","runtimes_o.txt"   ,sep=""))

#
###
