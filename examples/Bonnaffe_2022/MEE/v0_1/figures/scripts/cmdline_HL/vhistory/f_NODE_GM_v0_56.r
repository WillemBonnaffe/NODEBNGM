###############
## f_NODE_GM ##
###############

## goal: 
# - perfom residual minimisation training on the time series data
# - analyse the temporal change in variables with NODEs

## method:
# 1. interpolate the time series with sin ANN functions
# 2. estimate linear and non-linear coupling between the time series

## versions:
## 17-12-2020 - created version 0.0
## 18-12-2020 - created version 0.1
## 18-12-2020 - casted the model in a bayesian framework
## 18-12-2020 - created version 0.2
## 18-12-2020 - added a predictive model, before it was individual error models for each state variables
## 07-01-2021 - created version 0.3
## 07-01-2021 - checked code for potential errors 
## 07-01-2021 - created version 0.4
## 07-01-2021 - replaced polynomial process model by ANN process model 
## 08-01-2021 - created version 0.5
## 08-01-2021 - decoupled the interpolation from the fitting of the process model
## 14-01-2021 - implemented visualisation of the interpolation 
## 15-01-2021 - created version 0.6
## 15-01-2021 - implemented a bifurcation on the level of complexity considered
## 15-01-2021 - removed DEMCO as now obsolete because using gradients
## 15-01-2021 - computed r2 for each time series
## 15-01-2021 - created version 0.7
## 15-01-2021 - implemented computation of identifiability by ensembling
## 25-01-2021 - implemented computation of BIC of each model to formalise identification of complexity level 
## 30-01-2021 - implemented computation of the minimum acceptable prior complexity
## 03-02-2021 - created version 0.8
## 03-02-2021 - made sure that code was applicable to other time series by applying it to the tri-trophic system (Algee, Flagellate, Rotifer)
## 05-02-2021 - created version 0.9
## 05-02-2021 - implemented the bayesian regularisation scheme developped by Cawley and Talbot
## 14-02-2021 - created version 0.10
## 14-02-2021 - debugged code, switched back to bayesian regularisation with prior variance
## 15-02-2021 - implemented selection on complexity by minimising variance around estimates of effects through visual assessments
## 15-02-2021 - simplified code
## 22-02-2021 - created version 0.11
## 22-02-2021 - RE-implemented the Cawley approach
## 22-02-2021 - found a stable combination of regularisation parameters to fit the finch series with the standard normal approach
## 23-02-2021 - created version 0.12
## 23-02-2021 - tried to make the Cawley approach work 
## 24-02-2021 - created version 0.13
## 24-02-2021 - tried again the maximum constraint approach
## 24-02-2021 - created version 0.14
## 24-02-2021 - tested MCMC sampling instead to find the expected interpolation 
##            => it is not better, the model is still not fitted after 1,000,000 iteration
## 24-02-2021 - created version 0.15
##            - removed MCMC sampling
## 24-02-2021 - created version 0.16
##            - settled for the standard normal approach to regularisation (results will require a sensitivity analysis on prior sd)
##            - compared results coming from raw differences
## 24-02-2021 - created version 0.16
##            - settled for the standard normal approach to regularisation (results will require a sensitivity analysis on prior sd)
##            - compared results coming from raw differences
## 05-03-2021 - created version 0.17
##            - applied to the rotifer time series again
## 09-03-2021 - created version 0.18
##            - cleaned the code
##            - made code applicable to any time series
## 15-03-2021 - created version 0.19
##            - made figures for paper
## 17-03-2021 - polished figures
## 29-03-2021 - created version 0.20
##            - separated the fitting and the analysis/plotting
## 24-05-2021 - created version 0.21
##            - wrapped the approach in a single modular function
## 26-05-2021 - created version 0.22
##            - correct the calculation of the contribution matrix to no give weight to a variable when all contributions are small
## 27-05-2021 - created version 0.23
##            - implemented a cutoff of effects of the variables, if they are not significant they are set to be 0 
##            - changed the computation of contributions to be the sum of absolute values
## 27-05-2021 - created version 0.24
##            - compute significance of sum of squarred contributions vs a chi-square distribution
##            - moved the function .plot.DIN out of the function fit.NODE
## 01-06-2021 - created version 0.25 
##            - version now returns quantiles of mean effects and total contributions
## 09-06-2021 - created version 0.26
##            - polished figures 
## 14-06-2021 - created version 0.27
##            - simplified figures (e.g. removed the display of the time series and of net contributions in the summary figure)
##            - fit.NODE function now returns the ensembles for the quantity of interest (interpolation, effects, contributions, relative contributions)
## 18-11-2021 - created version 0_28
##            - renamed file f_NODE_GM_v0_28.r
##            - re-organised and cleaned code
## 22-11-2021 - cleaned code
## 22-11-2021 - created v0_29
##            - updated the main function
## 23-11-2021 - updated comments in code
## 23-11-2021 - created v0_30
##            - implemented scaling of process model
## 29-11-2021 - created v0_31
##            - cleaned code further
##            - removed marginal likelihood optimisation for simplicity 
##            - removed current global fitting function
##            - implemented marginal likelihood for interpolation as more stable
## 30-11-2021 - created v0_32
##            - cleaned code
##            - removed functions that are not essential to fitting NODE with GM
##            - implemented main function
## 18-02-2022 - created v0_33
##            - separated the fitting and the analysis of the NODE
## 23-02-2022 - created v0_34
##            - added function to sample NODE model
##            - cleaned and added code annotations
## 24-02-2022 - created v0_35
##            - cleaned code by removing dependency on global variables 
## 28-02-2022 - created v0_36
##            - changed yaxis in graph
## 03-03-2022 - created v0_37
##            - visualisation of results on natural scaled (not standardised)
## 03-03-2022 - created v0_38
##            - adapted code to work with variables with positive support (with no need for log)
## 14-03-2022 - created v0_39
##            - added functions to perform DFT
##            - added utility functions (for heatmaps etc.)
## 14-03-2022 - created v0_40
##            - updated main functions (e.g. fit.model_o)
## 15-03-2022 - created v0_41
##            - updated process model so that it can handle single covariates 
## 16-03-2022 - updated plot functions
##            - added main functions 
## 17-03-2022 - create v0_42
##            - added function to calculate relative contribution
##            - added function diamond plot
## 18-03-2022 - created v0_43
##            - improved plot functions
##            - created v0_44
##            - simplified and cleaned code where possible 
##            - instead of embedding scaling in ensemble return as metadata of the ensemble (along with time steps)
## 21-03-2022 - created v0_45
##            - simplified and cleaned code where possible 
## 21-03-2022 - created v0_46
##            - implemented polynomial perceptron
## 23-03-2022 - created v0_47
##            - cleaned code for perceptron
##            - re-implemented DEMC
## 25-03-2022 - created v0_48
##            - moved DEMC to annex function file
## 27-03-2022 - added legend in plot function
## 28-03-2022 - created v0_49
##            - added functions to calculate performance (likelihood, r-squared)
##            - created v0_50
##            - implemented cross validation 
##            - split the standardisation and fit
##            - created v0_51
## 08-04-2022 - created v0_52
##            - removed master function that are not required anymore
## 10-04-2022 - created v0_53
##            - moved initial parameter vector outside of fit function
##              => fix number of param issue in fit.model_p
##            - allow to switch activation functions
## 14-04-2022 - created v0_54
##            - functions predict of observation model also return ensemble predictions
## 25-04-2022 - created v0_55
##            - simplified the code further by removing plotting functions etc.
##            - combined bayesian model functions from observation and process model
## 25-05-2022 - created v0_56
##            - implemented weights for linear and nonlinear network

## next:
##            - make sure description of functions are accurate

###################
## FUNCTIONS SLP ##
###################

## activation functions
lin       = function(x) x
ddx.lin   = function(x) 1
pol2      = function(x) x^2
ddx.pol2  = function(x) 2*x
sigmo     = function(x) 1/(1+exp(-x)) 
ddx.sigmo = function(x) sigmo(x) * (1 - sigmo(x))
expo     = function(x) exp(x)
ddu.expo = function(x) exp(x)
# f_sigma     = function(x) exp(-x)
# ddu.f_sigma = function(x) -exp(-x)
# f_sigma     = function(x) (exp(2*x)-1)/(exp(2*x)+1) 
# ddu.f_sigma = function(x) (2*exp(2*x)-1)/(exp(2*x)+1) + (exp(2*x)-1)*(-1)*(2*exp(2*x)+1)/((exp(2*x)+1)^2)
# f_sigma     = function(x) sin(2*pi*x) 
# ddu.f_sigma = function(x) 2*pi*cos(2*pi*x) 
# f_sigma     = function(x) x^3
# ddu.f_sigma = function(x) 3*x^2
# f_sigma     = function(x) (x>0)*x
# ddu.f_sigma = function(x) (x>0)
# f_sigma     = function(x) exp(-x^2)
# ddu.f_sigma = function(x) -2*x*exp(-x^2)

## SLP ##
## goal: single layer perceptron returning single output with input x and parameter vector Omega 
# x       - vector - input variables
# Omega   - vector - parameters
# f_sigma - func - activation function
SLP = function(x,Omega,f_sigma) 
{	
	Omega = matrix(Omega,ncol=2 + length(x))
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
}

## ddOmega.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each parameter 
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddOmega.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
	Omega     = matrix(Omega,ncol=2 + length(x))
	x         = t(x)
	Omega_1   = Omega[,1]
	Omega_2   = Omega[,2]
	Omega_3   = Omega[,-c(1:2)]
	ddOmega_1 = f_sigma(Omega_2 + Omega_3%*%t(x))
	ddOmega_2 = Omega_1 * ddu.f_sigma(Omega_2 + Omega_3%*%t(x))
	ddOmega_3 = Omega_1%*%x * as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x)))
	return(c(ddOmega_1,ddOmega_2,ddOmega_3))
}

## ddx.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each input variable
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddx.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	x         = t(x)
	Omega_1   = Omega[,1]
	Omega_2   = Omega[,2]
	Omega_3   = Omega[,-c(1:2)]
	ddx = Omega_1%*%(Omega_3*as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x))))
	return(ddx)
}

#
###

##############################
## FUNCTIONS BAYESIAN MODEL ##
##############################

## goal: functions to define a simple Bayesian model with Gaussian error structure

## r2_p ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
r2 = function(X,Y,f,Omega)
{
  res      = Y - f(X,Omega) 
  r2       = 1 - sd(res^2)/sd(Y^2)
  return(r2)
}

## logLik ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
logLik = function(X,Y,f,Omega,sd_1)
{
  res      = Y - f(X,Omega)
  logLik   = - sum((res^2)/(sd_1^2))
  return(logLik)
}

## logPrior ##
## goal: compute log prior density of the process model 
# Omega - vector - parameters
# sd_2  - float  - standard deviation of prior
logPrior = function(Omega,sd_2)
{
  logPrior = - sum((Omega^2)/(sd_2^2))
  return(logPrior)
}

## logPost ##
## goal: log posterior distribution with normal error 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
logPost = function(X,Y,f,Omega,sd_1,sd_2)
{
	res      = Y - f(X,Omega)
	logLik   = - sum(  (res^2)/(sd_1^2))
	logPrior = - sum((Omega^2)/(sd_2^2))
	logPost  = logLik + logPrior
	return(logPost)
}

## ddOmega.logPost ##
## goal: compute the derivate of the log posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
ddOmega.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
	res              = Y - f(X,Omega)
	ddOmega.res      =   - df(X,Omega)
	ddOmega.logLik   = - 2 * ddOmega.res%*%res/(sd_1^2)
	ddOmega.logPrior = - 2 * Omega/(sd_2^2)
	ddOmega.logPost  = ddOmega.logLik + ddOmega.logPrior
	return(ddOmega.logPost)
}

## argmax.logPost ##
## goal: compute parameter vector that maximises log posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
argmax.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
	error_     = function(x) -logPost(X,Y,f,x,sd_1,sd_2)
	graderror_ = function(x) -ddOmega.logPost(X,Y,f,df,x,sd_1,sd_2)
	Omega      = optim(par    = Omega,
					   fn     = error_,
					   gr     = graderror_,
					   method = "BFGS"#,
					   # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
					   )$par
	return(Omega)
}

## logMar ##
## goal: compute the log marginal posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
logMar = function(X,Y,f,Omega)
{
	res      = Y - f(X,Omega)
	logLik   = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
	logPrior = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
	logMar   = logLik + 0.1*logPrior
	return(logMar)
}

## ddOmega.logMar ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
ddOmega.logMar = function(X,Y,f,df,Omega)
{
	res              = Y - f(X,Omega)
	ddOmega.res      =   - df(X,Omega)
	ddOmega.logLik   = - 0.5 * length(Y)     * 1/(0.5 * sum(res^2)   + 1) * 0.5 * ddOmega.res%*%res
	ddOmega.logPrior = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
	ddOmega.logMar   = ddOmega.logLik + 0.1*ddOmega.logPrior
	return(ddOmega.logMar)
}

## argmax.logMar ##
## goal: compute parameter vector that maximises the log marginal density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
argmax.logMar = function(X,Y,f,df,Omega)
{
	error_     = function(x) -logMar(X,Y,f,x)
	graderror_ = function(x) -ddOmega.logMar(X,Y,f,df,x)
	Omega      = optim(par    = Omega,
			           fn     = error_,
			           gr     = graderror_,
			           method = "BFGS"# ,
			           # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
			           )$par
	return(Omega)
}

#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## goal: functions for the observation model 

## next:
#       - move to loadModel_o
#       - find a way to do without predict

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

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## next:
#       - move to loadModel_p        
#       - find a way to do without predict

## activation function of process model
f_sigma_p     = expo     # sigmo     # lin     #  
ddx.f_sigma_p = ddu.expo # ddx.sigmo # ddx.lin #  
# f_sigma_p     = sigmo     # lin     #  
# ddx.f_sigma_p = ddx.sigmo # ddx.lin #  

# ## f_p ##
# ## goal: compute predicted response variable of process model defined given vector of explanatory variable at a given time step 
# # x       - vector - input variables
# # Omega   - vector - parameters
# f_p = function(x,Omega) 
# {	
# 	return(SLP(x,Omega,f_sigma_p))
# }
# 
# ## ddOmega.f_p ##
# ## goal: compute derivative vector of the process model wtr to each parameter at a given time step
# # x      - vector - input variables
# # Omega  - vector - parameters 
# ddOmega.f_p = function(x,Omega)
# {
# 	return(ddOmega.SLP(x,Omega,f_sigma_p,ddx.f_sigma_p))
# }
# 
# ## ddx.f_p ##
# ## goal: compute derivative vector of the process model wtr to each state variable at a given time step
# # x      - vector - input variables
# # Omega  - vector - parameters
# ddx.f_p = function(x,Omega)
# {
# 	return(ddx.SLP(x,Omega,f_sigma_p,ddx.f_sigma_p))
# }

w = 0.2

## f_p ##
## goal: compute predicted response variable of process model defined given vector of explanatory variable at a given time step 
# x       - vector - input variables
# Omega   - vector - parameters
f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2)
	return((1-w)*SLP(x,Omega[,1],lin) + (w)*SLP(x,Omega[,2],f_sigma_p))
}

## ddOmega.f_p ##
## goal: compute derivative vector of the process model wtr to each parameter at a given time step
# x      - vector - input variables
# Omega  - vector - parameters 
ddOmega.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return(c((1-w)*ddOmega.SLP(x,Omega[,1],lin,ddx.lin),(w)*ddOmega.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p)))
}

## ddx.f_p ##
## goal: compute derivative vector of the process model wtr to each state variable at a given time step
# x      - vector - input variables
# Omega  - vector - parameters
ddx.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2)
	return((1-w)*ddx.SLP(x,Omega[,1],lin,ddx.lin) + (w)*ddx.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p))
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
