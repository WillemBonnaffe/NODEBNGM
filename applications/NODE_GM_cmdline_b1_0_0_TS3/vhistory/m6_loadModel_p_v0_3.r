######################
## m6_loadModel_p.r ##
######################

## goal: load process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0
##            - allow modularity in process model (create new f_p function in loadModel script)
##            - fix number of parameters issue
## 14-04-2022 - created v0_1
##            - fitting of process model is done on a different ensemble member of Yhat_o at every iteration
## 15-04-2022 - simplified code
## 26-05-2022 - created v0_2
##            - moved functions of observation model to this script
##            - implemented weight parameter controlling for degree of nonlinearity
## 30-05-2022 - created v0_3
##            - implemented diffirential regularisation parameters
##              this is to allow reducing constraints on linear model
##              and increasing constraints on specific variables
##            - removed "slider"/weight parameter

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("f_NODE_GM_v0_58.r")
source("m5_loadData_p.r")

## parameters process model
K_p   = 10
W_p   = c(10,10,10)
w     = 0.4

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## activation function of process model
f_sigma_p     = expo     
ddx.f_sigma_p = ddu.expo 
# f_sigma_p     = sigmo     
# ddx.f_sigma_p = ddx.sigmo 

## f_p ##
## goal: compute predicted response variable of process model defined given vector of explanatory variable at a given time step 
# x       - vector - input variables
# Omega   - vector - parameters
f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2)
	return(SLP(x,Omega[,1],lin) + w*SLP(x,Omega[,2],f_sigma_p))
}

## ddOmega.f_p ##
## goal: compute derivative vector of the process model wtr to each parameter at a given time step
# x      - vector - input variables
# Omega  - vector - parameters 
ddOmega.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return(c(ddOmega.SLP(x,Omega[,1],lin,ddx.lin),w*ddOmega.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p)))
}

## ddx.f_p ##
## goal: compute derivative vector of the process model wtr to each state variable at a given time step
# x      - vector - input variables
# Omega  - vector - parameters
ddx.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return(ddx.SLP(x,Omega[,1],lin,ddx.lin) + w*ddx.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p))
}
 
## *.eval ##
## goal: compute process model functions across several time steps
# X      - matrix - matrix of variables across all time steps
# ddt.X  - matrix - matrix of temporal derivatives of variables
# Omega  - vector - vector of parameters of the network 
f_p.eval           = function(X,Omega) apply(t(X),2,function(x) f_p(x,Omega))
ddx.f_p.eval       = function(X,Omega) apply(t(X),2,function(x) ddx.f_p(x,Omega))
ddOmega.f_p.eval   = function(X,Omega) apply(t(X),2,function(x) ddOmega.f_p(x,Omega))

#
###
