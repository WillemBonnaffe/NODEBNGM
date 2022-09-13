######################
## m5_loadModel_p.r ##
######################

## goal: load process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0
##            - allow modularity in process model (create new f_p function in loadModel script)
##            - fix number of parameters issue

## next:
##            - evaluate predictions for all Omega_iku
##            - parallelise code

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## imports
source("f_NODE_GM_v0_53.r")
source("m1_loadData.r")

## load results observation model
load(paste(pathToOut,"/","Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))

## convert list to table
Yhat_o_     = NULL
ddt.Yhat_o_ = NULL
for(i in 1:length(Yhat_o))
{   
    Yhat_o_     = cbind(Yhat_o_,       Yhat_o[[i]][["MaP.Yhat_o"]])
    ddt.Yhat_o_ = cbind(ddt.Yhat_o_,ddt.Yhat_o[[i]][["MaP.ddt.Yhat_o"]])
}

## prepare data
X_p     = Yhat_o_
ddt.X_p = ddt.Yhat_o_ 
Y_p     = 1/Yhat_o_ * ddt.Yhat_o_
colnames(Y_p) = colnames(TS)[-1]

## parameters process model
K_p   = 10                  # number of ensemble elements
W_p   = 10                  # number of nodes in process model
N_p   = 2*W_p*(2+ncol(X_p))
sd1_p = 0.1                 # sd in likelihood of residuals of process model
sd2_p = seq(0.05,0.5,0.05)  
folds = list(c(0,1/3),c(1/3,2/3),c(2/3,3/3)) # list(c(0,2/3),c(0.5/3,2.5/3),c(1/3,3/3)) # list(c(0,2/4),c(1/4,3/4),c(2/4,4/4))

#
###
