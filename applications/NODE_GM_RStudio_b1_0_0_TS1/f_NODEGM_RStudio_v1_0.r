##############
## f_NODEGM ##
##############

## goal: 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## method:
## 1. interpolate the time series with sin ANN functions
## 2. estimate linear and non-linear coupling between the time series

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
## 26-05-2022 - created v0_57
##            - moved model o and model p functions to external scripts
##            - cleaned code 
## 30-05-2022 - created v0_58
##            - implemented compensation parameter "c" to counteract the increase in network size in loss function of marginal posterior
##            - removed slider/weight parameter as not necessary anymore given that approach is now to vary size of network
## 09-06-2022 - created v0_59
##            - introduced functions for fitting/training and plotting to the main function file
##              this is to make it easier to deploy the NODEGM approach
## 09-06-2022 - renamed f_NODEGM_v0_59 to f_NODEGM_compact_v1_0
##            - this is the version for release


###################
## FUNCTIONS SLP ##
###################

## activation functions ##
lin       = function(x) x
ddx.lin   = function(x) 1
pol2      = function(x) x^2
ddx.pol2  = function(x) 2*x
sigmo     = function(x) 1/(1+exp(-x)) 
ddx.sigmo = function(x) sigmo(x) * (1 - sigmo(x))
expo      = function(x) exp(x)
ddx.expo  = function(x) exp(x)
relu      = function(x) (x>0)*x
ddx.relu  = function(x) (x>0)
tanh      = function(x) (exp(2*x)-1)/(exp(2*x)+1) 
ddx.tanh  = function(x) (2*exp(2*x))/(exp(2*x)+1) + (exp(2*x)-1)*(-1)*(2*exp(2*x))/((exp(2*x)+1)^2)

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

## logMarLik ##
## goal: compute the log marginal likelihood 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarLik = function(X,Y,f,Omega)
{
	res      = Y - f(X,Omega)
	logMarLik   = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
	return(logMarLik)
}

## logMarPri ##
## goal: compute the log marginal prior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarPri = function(Omega)
{
	logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1) 
    # Ntmp = length(Omega)
	# logMarPri = - 0.5 * length(Omega[(Ntmp/2+1):Ntmp]) * log(0.5 * sum(Omega[(Ntmp/2+1):Ntmp]^2) + 1)
	return(logMarPri)
}

## logMarPost ##
## goal: compute the log marginal posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarPost = function(X,Y,f,Omega,c=1)
{
	res       = Y - f(X,Omega)
	logMarLik = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
	logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
    # Ntmp = length(Omega)
	# logMarPri = - 0.5 * length(Omega[(Ntmp/2+1):Ntmp]) * log(0.5 * sum(Omega[(Ntmp/2+1):Ntmp]^2) + 1)
	logMarPos = logMarLik + c*logMarPri
	return(logMarPos)
}

## ddOmega.logMarPost ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - parameters
ddOmega.logMarPost = function(X,Y,f,df,Omega,c=1)
{
res                = Y - f(X,Omega)
	ddOmega.res        =   - df(X,Omega)
	ddOmega.logMarLik  = - 0.5 * length(Y)     * 1/(0.5 * sum(res^2)   + 1) * 0.5 * ddOmega.res%*%res
	ddOmega.logMarPri  = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
    # Ntmp = length(Omega)
	# ddOmega.logMarPri  = - 0.5 * length(Omega[(Ntmp/2+1):Ntmp]) * 1/(0.5 * sum(Omega[(Ntmp/2+1):Ntmp]^2) + 1) * Omega[(Ntmp/2+1):Ntmp]
	ddOmega.logMarPos  = ddOmega.logMarLik + c*ddOmega.logMarPri ## divide by number of neurones in the network
	return(ddOmega.logMarPos)
}

## argmax.logMarPost ##
## goal: compute parameter vector that maximises the log marginal density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - fn     - function to predict response
# df    - fn     - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
argmax.logMarPost = function(X,Y,f,df,Omega,c=1)
{
	error_     = function(x) -logMarPost(X,Y,f,x,c)
	graderror_ = function(x) -ddOmega.logMarPost(X,Y,f,df,x,c)
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

## loadData_o ##
## goal: prepare data for training of observation model
loadData_o = function(TS,alpha_i)
{
    ## data specs
    N       = ncol(TS) - 1
    n       = nrow(TS)
        
    ## predictive and response variable
    t   = TS[,1] 
    Y   = TS[,-1]
    nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)
    
    ## standardise time steps
    t_0 = min(t)
    t_f = max(t)
    dt  = diff(t[1:2])
    t_  = 1:length(t)
    nt_ = seq(min(t_),max(t_),(t_[2]-t_[1])/alpha_i)
    
    ## standardise data
    Y_     = log(Y)
    mean_y = apply(Y_,2,mean)
    sd_y   = apply(Y_,2,sd)
    Y_     = t((t(Y_)-mean_y)/sd_y)

    ## terminate
    return(list("t"      = t,
                "t_"     = t_,
                "nt"     = nt,
                "nt_"    = nt_,
                "Y"      = Y,
                "Y_"     = Y_,
                "mean_y" = mean_y,
                "sd_y"   = sd_y,
                "dt"     = dt,
                "N"      = N
                ))
}

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

## train_o ##
trainModel_o = function(TS,alpha_i,W_o,K_o,rho=1)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)

    ## interpolate each column in the time series
    Omega_o    = list()
    Yhat_o     = list()
    ddt.Yhat_o = list()
    logPost_o  = list()
    for (i in 1:N)
    {
     	## iterator
     	message(paste("fitting: ",i,"/",N,sep=""))
    
    	## get ensemble
    	Omega_o_i = NULL
        logPost_o_i = NULL
    	for(k in 1:K_o)
    	{   
            check = F
            while(check == F)
            {
    		    ## fit
    	        Omega_0   = rnorm(W_o[i]*3,0,0.001)
    	        Omega_f   = argmax.logMarPost(t_,Y_[,i],f_o.eval,ddOmega.f_o.eval,Omega_0,1/W_o[i])
    
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
                    logPost_o_i = c(logPost_o_i,logMarLik(t_,Y_[,i],f_o.eval,Omega_f))
                }
            }
    	}
    
        ## sort by best models
        s = rev(order(logPost_o_i))[1:round(rho*K_o)]
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
    return(list("Yhat_o"     = Yhat_o,
                "ddt.Yhat_o" = ddt.Yhat_o,
                "Omega_o"    = Omega_o))
}

## plot_o ## 
plotModel_o = function(TS,alpha_i=1,Yhat_o,ddt.Yhat_o)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)

    col   = rev(rainbow(N,start=0.6,end=0.9))
    xlab  = c("","Time")
    ylab  = c("Y(t)","dY/dt(t)")
    index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
    par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
    layout(cbind(c(1,2),c(3,4),c(5,6))[,1:min(3,N)])
    for (i in 1:N)
    {
    	## predictions 
        E.Yhat_o       = apply(Yhat_o[[i]],2,mean)
        q05.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.05)
        q95.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.95)
        E.ddt.Yhat_o   = apply(ddt.Yhat_o[[i]],2,mean)
        q05.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.05)
        q95.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.95)
        #
    	## training data
    	t = TS[,1] 
    	Y = TS[,-1][,i] 
    	nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)
    	#
    	## visualise interpolated response
    	plot(t,rep(0,length(t)),ylim=c(min(Y)-0.2,max(Y))+0.2*c(-min(Y),max(Y)),type="l",lty=3,xlab=xlab[1],ylab=if(i == 1)ylab[1]else"")
    	polygon(c(nt,rev(nt)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
    	lines(nt,E.Yhat_o,col=adjustcolor(col[i],0.75),lwd=2)
    	points(t,Y,pch=16,col=adjustcolor("black",0.75)) 
    	if(!is.null(index)) legend("topright",legend=index[1+(i-1)*2],bty="n",cex=1.5)
        legend("bottom",legend=c("G. truth","Estimate"),lty=c(2,1),col=c(adjustcolor("black",0.75),col[i]),lwd=2,bty="n",horiz=T)
    	#
    	## visualise temporal derivative
    	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*2,type="l",lty=3,xlab=xlab[2],ylab=if(i == 1)ylab[2]else"")
    	polygon(c(nt,rev(nt)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
    	lines(nt,E.ddt.Yhat_o,col=adjustcolor(col[i],0.75),lwd=2)
    	if(!is.null(index)) legend("topright",legend=index[2+(i-1)*2],bty="n",cex=1.5)
        legend("bottom",legend=c("G. truth","Estimate"),lty=c(2,1),col=c(adjustcolor("black",0.75),col[i]),lwd=2,bty="n",horiz=T)
    }
}

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## loadData_p ##
loadData_p = function(Yhat_o,ddt.Yhat_o)
{
    ## data specs
    N   = length(Yhat_o)
    n   = ncol(Yhat_o[[1]])
    K_o = nrow(Yhat_o[[1]])

    ## predictive variable
    E.X_p     = NULL
    E.ddt.X_p = NULL
    E.Y_p     = NULL
    X_p       = NULL
    ddt.X_p   = NULL
    Y_p       = NULL
    for(j in 1:N)
    {
        # ## select random ensemble element (for stochastic training)
        # v       = sample(1:K_o,1)
        # X_p     = cbind(X_p,        Yhat_o[[j]][v,])
        # ddt.X_p = cbind(ddt.X_p,ddt.Yhat_o[[j]][v,])
        # Y_p     = cbind(Y_p,1/Yhat_o[[j]][v,]*ddt.Yhat_o[[j]][v,])
    
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

    ## terminate
    return(list("X"      = X,
                "ddt.X"  = ddt.X,
                "X_"     = X_,
                "Y"      = Y,
                "Y_"     = Y_,
                "mean_x" = mean_x,
                "sd_x"   = sd_x,
                "mean_y" = mean_y,
                "sd_y"   = sd_y
                ))
}

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

## trainModel_p ##
trainModel_p = function(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p)
{
    ## for each variable
    Omega_p    = list()
    Yhat_p     = list()
    ddx.Yhat_p = list()
    Geber_p    = list()
    for(i in 1:N)
    {
        ## iterator
        message(paste("fitting: ",i,"/",N,sep=""))
    
        ## fit model
        Omega_p_i = NULL
        for(k in 1:K_p)
        {
    
            ## dataloader
            attach(loadData_p(Yhat_o,ddt.Yhat_o),warn.conflicts=F)
    
            ## fit
            Omega_0      = rnorm(N_p[i],0,0.001)
            Yhat         = function(X,Omega) f_p.eval(X,Omega)
            ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega)
            # Omega_f      = argmax.logMarPost(X_,Y_[,i],Yhat,ddOmega.Yhat,Omega_0,1/W_p[i])
            Omega_f      = argmax.logPost(X_,Y_[,i],Yhat,ddOmega.Yhat,Omega_0,sd1_p,sd2_p[[i]])
            Omega_p_i    = rbind(Omega_p_i,Omega_f)
    
            ## update
            logPost_0    = logLik(X_,Y_[,i],Yhat,Omega_0,sd1_p)
            logPost_f    = logLik(X_,Y_[,i],Yhat,Omega_f,sd1_p)
            # logPost_0    = logMarLik(X_,Y_[,i],Yhat,Omega_0)
            # logPost_f    = logMarLik(X_,Y_[,i],Yhat,Omega_f)
            message(paste(k,"/",K_p,"\t",
                    format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
                    format(round(logPost_f,2),nsmall=2),sep=""))
        }
    
        ## store
        Omega_p[[i]] = Omega_p_i
    
        ## predictions
        Yhat_p[[i]]     = t(apply(Omega_p[[i]],1,function(x) mean_y[i] + sd_y[i] * f_p.eval(X_,x)))
        ddx.Yhat_p[[i]] = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * ddx.f_p.eval(X_,x)))
        Geber_p[[i]]    = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * t(ddt.X) * ddx.f_p.eval(X_,x)))
    }

    ## terminate
    return(list("Yhat_p"     = Yhat_p,
                "ddx.Yhat_p" = ddx.Yhat_p,
                "Geber_p"    = Geber_p,
                "Omega_p"    = Omega_p
                ))
}

## plotModel_p ##
plotModel_p = function(TS,alpha_i=1,Yhat_p,ddx.Yhat_p,Geber_p)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)
    attach(loadData_p(Yhat_o,ddt.Yhat_o),warn.conflicts=F)
    
    col    = rev(rainbow(N,start=0.6,end=0.9))
    xlab   = c("","","Time")
    ylab   = c("P.c. growth rate","Effects","Contributions")
    index  = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
    legend = paste(colnames(TS)[-1])
    par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
    layout(cbind(c(1,2,3),c(4,5,6),c(7,8,9))[,1:min(3,N)])
    for(i in 1:N)
    {
        ## predictions
        E.Yhat_p       = apply(Yhat_p[[i]],2,mean)
        q05.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.05)
        q95.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.95)
        E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p[[i]],2,mean),ncol=nrow(X)))
        q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
        q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
        E.Geber_p      = t(matrix(apply(Geber_p[[i]],2,mean),ncol=nrow(X)))
        q05.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
        q95.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
        #
        ## dynamics
        x = nt
        y = Y[,i]
        plot(x,rep(0,length(x)),ylim=c(-1,1)*max(abs(y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=if(i==1)ylab[1]else"")
        points(x,y,pch=16,col=adjustcolor("black",0.75))
        polygon(c(nt,rev(nt)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col[i],alpha=0.2),border=NA)
        lines(nt,E.Yhat_p,col=adjustcolor(col[i],alpha=0.75),lwd=2)
        if(!is.null(index)) legend("topright",legend=index[1+(i-1)*(3)],bty="n",cex=1.5)
        #
        ## effects
        plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=if(i==1)ylab[2]else"")
        for(j in 1:N) lines(nt,E.ddx.Yhat_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
        for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
        if(!is.null(index))  legend("topright",legend=index[2+(i-1)*(3)],bty="n",cex=1.5)
        if(!is.null(legend)) legend("bottom" ,legend=legend,lty=1,col=adjustcolor(col,alpha=0.75),bty="n",horiz=T,lwd=2)
        #
        ## Geber
        plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=if(i==1)ylab[3]else"")
        for(j in 1:N) lines(nt, E.Geber_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
        for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
        if(!is.null(index))  legend("topright",legend=index[3+(i-1)*(3)],bty="n",cex=1.5)
        if(!is.null(legend)) legend("bottom" ,legend=legend,lty=1,col=adjustcolor(col,alpha=0.75),bty="n",horiz=T,lwd=2)
    }
}

## crossVal_p ##
crossVal_p = function(TS,alpha_i=1,Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,folds,crossValParVect)
{
    ## data specs
    N   = length(Yhat_o)
    n   = ncol(Yhat_o[[1]])

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
                    attach(loadData_o(TS,alpha_i),warn.conflicts=F)
    
                    ## split training/test
                    s   = round(folds[[u]][1]*n+1):round(folds[[u]][2]*n)
                    X_l = X_[-s,]
                    Y_l = Y_[-s,]
                    X_t = X_[s,]
                    Y_t = Y_[s,]
    
                    ## TO MODULARISE ##
                    sd2_p[[i]][(N_p[i]/2):N_p[i]] = crossValParVect[k] 
    
                    ## fit
                    Omega_0      = rnorm(N_p[i],0,0.001)
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
                    #       format(round(logMarPost_0,2),nsmall=2),"\t","-->","\t",
                    #       format(round(logMarPost_f,2),nsmall=2),sep=""))
    
                    ## model performance
                    logLik_l = logLik(X_l,Y_l[,i],Yhat,Omega_f,sd1_p)
                    logLik_t = logLik(X_t,Y_t[,i],Yhat,Omega_f,sd1_p)
                    # logMarLik_l = logMarLik(X_l,Y_l[,i],Yhat,Omega_f)
                    # logMarLik_t = logMarLik(X_t,Y_t[,i],Yhat,Omega_f)
                    # message(paste(m,"/",K_p,"\t",
                    #       format(round(logMarLik_l,2),nsmall=2),"\t","","\t",
                    #       format(round(logMarLik_t,2),nsmall=2),sep=""))
    
                    ## cross validation matrix
                    crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))
                    # crossVal_ik  = rbind(crossVal_ik,cbind(logMarLik_l,logMarLik_t))
    
                }
                Omega_p[[i]][[k]][[u]] = Omega_p_iku
            }
            message("")
    
            ## store
            E.crossVal_ik  = apply(crossVal_ik,2,mean)
            sd.crossVal_ik = apply(crossVal_ik,2,sd)
            crossVal_i     = rbind(crossVal_i,c(crossValParVect[k],E.crossVal_ik,sd.crossVal_ik))
            message(paste("logLik l vs t: ",
                    format(round(E.crossVal_ik[1],2),nsmall=2),"\t",
                    format(round(E.crossVal_ik[2],2),nsmall=2),sep=""))
    
        }
        message("\n")
    
        ## store
        crossVal_p[[i]] = crossVal_i
        colnames(crossVal_p[[i]]) = c("sd","logLik_l","logLik_t","sd.logLik_l","sd.logLik_t")
        # colnames(crossVal_p[[i]]) = c("w","logMarLik_l","logMarLik_t","sd.logMarLik_l","sd.logMarLik_t")
    }
    message("\n")
    return(crossVal_p)
}

## figure cross validation ##
plotCrossVal_p = function(resultsCrossVal)
{
    index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
    colVect = rainbow(2,start=0.6,end=0.9)
    par(mfrow=c(min(3,N),1),mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
    for(i in 1:N)
    {
        ## unpack
        crossVal = resultsCrossVal_p[[i]]
        #
        ## plot
        x    = crossVal[,"sd"]
        y    = crossVal[,c("logLik_l","logLik_t")]
        sd.y = crossVal[,c("sd.logLik_l","sd.logLik_t")]
        xlab = if(i == N) "Std. prior" else ""
        ylab = "Log likelihood"
        plot(x,rep(0,length(x)),ylim=c(min(y-sd.y)-0.5*abs(min(y-sd.y)),max(y+sd.y)+0.5*abs(max(y+sd.y))),type="l",lty=3,xlab=xlab,ylab=ylab)
        #
        ## training line
        col  = colVect[1]
        y    = crossVal[,"logLik_l"]
        sd.y = crossVal[,"sd.logLik_l"]
        lines(x,y,col=col)
        polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
        #
        ## testing line
        col  = colVect[2]
        y    = crossVal[,"logLik_t"]
        sd.y = crossVal[,"sd.logLik_t"]
        lines(x,y,col=col)
        polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
        #
        ## legend
        legend("bottom",legend=c("Training","Testing"),col=colVect,lty=1,bty="n",cex=1.5,horiz=T)
        legend("topright",legend=index[i],bty="n",cex=1.5)
    }
}

#
###
