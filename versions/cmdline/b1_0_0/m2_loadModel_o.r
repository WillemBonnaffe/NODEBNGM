#######################
## m2_loadModel_o.r ##
#######################

## goal: load observation model

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM.r")

## parameters observation model
K_o      = 10
W_o      = rep(100,N)
N_o      = W_o * 3
sd1_o    = 0.1
sd2_o    = rep(0.001,N) 
rho      = 1 

#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## goal: functions for the observation model 

## f_o ##
## goal: compute predicted values of response variable at time step t
# t     - float  - time step 
# Omega - vector - parameters 
f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o ##
## goal: compute time derivative of the predicted response t time step t
# t     - float  - time step 
# Omega - vector - parameters 
ddt.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega.f_o ##
## goal: compute derivative of the predicted response wtr to each network parameter
# t     - float  - time step
# Omega - vector - parameters 
ddOmega.f_o = function(t,Omega)
{	
	Omega      = matrix(Omega,ncol=3)
	dfdOmega_1 =                sin(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega_2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega_3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
	return(c(dfdOmega_1,dfdOmega_2,dfdOmega_3))
}

## *.eval ##
## goal: compute functions across multiple time steps
# t     - vector - time steps in arbitrary units
# Omega - vector - parameters 
f_o.eval         = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval     = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega.f_o(x,Omega))

#
###


