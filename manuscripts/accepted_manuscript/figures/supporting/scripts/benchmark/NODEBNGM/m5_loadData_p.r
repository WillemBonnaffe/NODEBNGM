#####################
## m5_loadData_p.r ##
#####################

## goal: load data for training process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 15-04-2022 - created v0_0
## 25-04-2022 - created v0_1
##            - loading and processing ensemble instead of pre-computed mean

##############
## INITIATE ##
##############

## goal: initiate data for training process model

## make out directory
pathToOut = "out"

## load results observation model
load(paste(pathToOut,"/","Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))

## data specs 
N   = length(Yhat_o)
n   = ncol(Yhat_o[[1]])
K_o = nrow(Yhat_o[[1]])

#
###

##################
## PREPARE DATA ##
##################

## predictive variable
E.X_p     = NULL
E.ddt.X_p = NULL
E.Y_p     = NULL
X_p       = NULL
ddt.X_p   = NULL
Y_p       = NULL
for(j in 1:N)
{   
	## select random ensemble element
    v       = sample(1:K_o,1)
    X_p     = cbind(X_p,        Yhat_o[[j]][v,])
    ddt.X_p = cbind(ddt.X_p,ddt.Yhat_o[[j]][v,])
    Y_p     = cbind(Y_p,1/Yhat_o[[j]][v,]*ddt.Yhat_o[[j]][v,])

	## mean of ensemble
    E.X_p     = cbind(E.X_p,     apply(Yhat_o[[j]],2,mean))
    E.ddt.X_p = cbind(E.ddt.X_p, apply(ddt.Yhat_o[[j]],2,mean))
    E.Y_p     = cbind(E.Y_p,     apply(1/Yhat_o[[j]]*ddt.Yhat_o[[j]],2,mean))
}   

## variables
X     = E.X_p
ddt.X = E.ddt.X_p
Y     = E.Y_p

## standardise predictive variables
X_     = X
mean_x = apply(E.X_p,2,mean)
sd_x   = apply(E.X_p,2,sd)
X_     = t((t(X_)-mean_x)/sd_x)

## standardise response variable
Y_     = Y
mean_y = apply(E.Y_p,2,mean)
sd_y   = apply(E.Y_p,2,sd)
Y_     = t((t(Y_)-mean_y)/sd_y)

#
###
