######################
## m6_loadModel_p.r ##
######################

## goal: load process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("f_NODE_GM.r")
source("m5_loadData_p.r")

## parameters process model
t_train = 1:(round(2/3*n))
t_test  = (round(2/3*n)+1):n
K_p   = 30
W_p   = rep(10,N)
N_p   = 2 * W_p * (2+N)
sd1_p = 0.1
sd2_p = list(0.05,0.05,0.05)

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## activation function of process model
# f_sigma_p     = pol2   
# ddx.f_sigma_p = ddx.pol2
f_sigma_p     = expo     
ddx.f_sigma_p = ddx.expo 
# f_sigma_p     = sigmo     
# ddx.f_sigma_p = ddx.sigmo 

## f_p ##
## goal: compute predicted response variable of process model defined given vector of explanatory variable at a given time step 
# x       - vector - input variables
# Omega   - vector - parameters
f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2)
	return(SLP(x,Omega[,1],lin) + SLP(x,Omega[,2],f_sigma_p))
}

## ddOmega.f_p ##
## goal: compute derivative vector of the process model wtr to each parameter at a given time step
# x      - vector - input variables
# Omega  - vector - parameters 
ddOmega.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return(c(ddOmega.SLP(x,Omega[,1],lin,ddx.lin),ddOmega.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p)))
}

## ddx.f_p ##
## goal: compute derivative vector of the process model wtr to each state variable at a given time step
# x      - vector - input variables
# Omega  - vector - parameters
ddx.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return(ddx.SLP(x,Omega[,1],lin,ddx.lin) + ddx.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p))
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
