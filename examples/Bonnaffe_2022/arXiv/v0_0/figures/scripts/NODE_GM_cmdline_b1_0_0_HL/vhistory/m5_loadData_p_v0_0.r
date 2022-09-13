#####################
## m5_loadData_p.r ##
#####################

## goal: load data for training process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 15-04-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate data for training process model

## make out directory
pathToOut = "out/RPS/"

## load results observation model
load(paste(pathToOut,"/","Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))

## data specs 
N   = length(Yhat_o)
n   = ncol(Yhat_o[[1]][["Yhat_o"]])
K_o = nrow(Yhat_o[[1]][["Yhat_o"]])

#
###

##################
## PREPARE DATA ##
##################

## predictive variable
X_p     = NULL
ddt.X_p = NULL
for(j in 1:N)
{   
    v       = sample(1:K_o,1)
    X_p     = cbind(X_p,        Yhat_o[[j]][["Yhat_o"]][v,])
    ddt.X_p = cbind(ddt.X_p,ddt.Yhat_o[[j]][["ddt.Yhat_o"]][v,])
    # X_p     = cbind(X_p,        Yhat_o[[j]][["E.Yhat_o"]])
    # ddt.X_p = cbind(ddt.X_p,ddt.Yhat_o[[j]][["E.ddt.Yhat_o"]])
}   

## response variable
Y_p     = 1/X_p * ddt.X_p
colnames(Y_p) = names(Yhat_o)

## variables
X     = X_p
ddt.X = ddt.X_p
Y     = Y_p

## standardise predictive variables
X_     = X
mean_x = apply(X_,2,mean)
sd_x   = apply(X_,2,sd)
X_     = t((t(X_)-mean_x)/sd_x)

## standardise response variable
Y_     = Y
mean_y = apply(Y_,2,mean)
sd_y   = apply(Y_,2,sd)
Y_     = t((t(Y_)-mean_y)/sd_y)

#
###
