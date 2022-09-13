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

## next:
##            - fix number of param issue in fit.model_p
##            - add log of response variable in process model
##            - make sure description of functions are accurate
##            - review name of functions
##            - allow for different values of priors for each time series
##            - fix number of parameters in model
##            - allow to switch activation functions

###################
## DFT FUNCTIONS ##
###################

## goal: functions to perform discrete fourier transforms of time series

## DFT ##
## goal: compute discrete Fourier transform of a signal f
# f  - vector - vector of values of signal
# x  - vector - vector of dependent variable (e.g. time)
# x_ - vector - vector of dependent variable values at which to interpolate the signal
# K  - int    - number of elements in the series
DFT = function(f,x,x_,K)
{
    dx = diff(x[1:2])
    L  = max(x)
    n  = length(x)
    A0 = sum(f*rep(1,length(x)))*dx*2/L
    fFS = A0/2
    for (k in 1:K)
    {
        Ak = sum(f*cos(2*pi*k*x/L))*dx*2/L
        Bk = sum(f*sin(2*pi*k*x/L))*dx*2/L
        fFS = fFS + Ak*cos(2*pi*k*x_/L) + Bk*sin(2*pi*k*x_/L)
    }
    return(fFS)
}

## dDFT ##
## goal: compute the derivative of discrete Fourier transform of a signal f wtr to dependent variable x
# f  - vector - vector of values of signal
# x  - vector - vector of dependent variable (e.g. time)
# x_ - vector - vector of dependent variable values at which to interpolate the signal
# K  - int    - number of elements in the series
dDFT = function(f,x,x_,K)
{
    dx = diff(x[1:2])
    L  = max(x)
    n  = length(x)
    A0 = 0
    dfFS = A0/2
    for (k in 1:K)
    {
        Ak = sum(f*cos(2*pi*k*x/L))*dx*2/L
        Bk = sum(f*sin(2*pi*k*x/L))*dx*2/L
        dfFS = dfFS - Ak*2*pi*k/L*sin(2*pi*k*x_/L) + Bk*2*pi*k/L*cos(2*pi*k*x_/L)
    }
    return(dfFS)
}

## fit.DFT ##
## goal: compute interpolation/extrpolation of a time series
# Y   - vector - vector of values of time series at time steps x
# x   - vector - vector of time steps
# x_  - vector - vector of interpolated/extrapolated time steps
# K   - int    - number of elements in the series
# log - bool   - whether to log the values of the response Y (e.g. if strictly positive)
fit.DFT = function(Y,x,x_,K,log=F)
{
	## standardise timesteps [min(t),max(t)] -> [1,N]
	dt = diff(x[1:2])
	t  = seq(1,length(x))
	nt = (x_-min(x))/(max(x)-min(x))*(length(x)-1) + 1

	## scale data
	if (log==T) Y_ = log(Y) else Y_ = Y
	mean_ = mean(Y_)
	sd_   = sd(Y_)
	Y_    = (Y_-mean_)/sd_
	
	## compute DFT
	Yhat     = DFT(Y_,t,nt,K)
	ddt.Yhat = dDFT(Y_,t,nt,K)
	
	## de-scale data
	if (log==T)
	{	
		Yhat     = exp(mean_ + sd_ * Yhat)
		ddt.Yhat = 1/dt * sd_ * Yhat * ddt.Yhat
	} else
	{
		Yhat     = mean_ + sd_ * Yhat
		ddt.Yhat = 1/dt * sd_ * ddt.Yhat
	}

	## terminate
	return(cbind(Yhat,ddt.Yhat))
}

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

## logPost_o ##
## goal: compute the log posterior of the observation model given observed response 
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
logPost_o = function(t,Y,Omega,sd_1,sd_2)
{
	Yhat     = f_o.eval(t,Omega)
	res      = Y - Yhat
	logLik   = - sum(  (res^2)/(sd_1^2))
	logPrior = - sum((Omega^2)/(sd_2^2))
	logPost  = logLik + logPrior
	return(logPost)
}

## ddOmega.logPost_o ##
## goal: compute the derivate of the log posterior density wtr to each parameter
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
ddOmega.logPost_o = function(t,Y,Omega,sd_1,sd_2)
{
	Yhat             = f_o.eval(t,Omega)
	res              = Y - Yhat
	ddOmega.res      = - ddOmega.f_o.eval(t,Omega)
	ddOmega.logLik   = - 2 * ddOmega.res%*%res/(sd_1^2)
	ddOmega.logPrior = - 2 * Omega/(sd_2^2)
	ddOmega.logPost  = ddOmega.logLik + ddOmega.logPrior
	return(ddOmega.logPost)
}

## argmax.logPost_o ##
## goal: compute parameter vector that maximises log posterior density
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - initial parameters 
argmax.logPost_o = function(t,Y,Omega,sd_1,sd_2)
{
	error_     = function(x) -logPost_o(t,Y,x,sd_1,sd_2)
	graderror_ = function(x) -ddOmega.logPost_o(t,Y,x,sd_1,sd_2)
	Omega      = optim(par    = Omega,
					   fn     = error_,
					   gr     = graderror_,
					   method = "BFGS"#,
					   # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
					   )$par
	return(Omega)
}

## logMar_o ##
## goal: compute the log marginal posterior density of observation model
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - parameters
logMar_o = function(t,Y,Omega)
{
	Yhat     = f_o.eval(t,Omega)
	res      = Y - Yhat
	logLik   = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
	logPrior = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
	logMar   = logLik + logPrior
	return(logMar)
}

## ddOmega.logMar_o ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - parameters
ddOmega.logMar_o = function(t,Y,Omega)
{
	Yhat             = f_o.eval(t,Omega)
	res              = Y - Yhat
	ddOmega.res      = - ddOmega.f_o.eval(t,Omega)
	ddOmega.logLik   = - 0.5 * length(Y)     * 1/(0.5 * sum(res^2)   + 1) * 0.5 * ddOmega.res%*%res
	ddOmega.logPrior = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
	ddOmega.logMar   = ddOmega.logLik + ddOmega.logPrior
	return(ddOmega.logMar)
}

## argmax.logMar_o
## goal: compute parameter vector that maximises the log marginal density
# t     - vector - time steps 
# Y     - vector - observations
# Omega - vector - initial parameters 
argmax.logMar_o = function(t,Y,Omega)
{
	error_     = function(x) -logMar_o(t,Y,x)
	graderror_ = function(x) -ddOmega.logMar_o(t,Y,x)
	Omega      = optim(par    = Omega,
			           fn     = error_,
			           gr     = graderror_,
			           method = "BFGS"# ,
			           # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
			           )$par
	return(Omega)
}
	
## fit.model_o ##
## goal: fit the observation model (i.e. interpolate/extrapolate the time series)
## out:  returns table of ensemble of parameter vectors 
# t      - vector - time steps
# Y      - vector - observations 
# W      - int    - number of hidden nodes in the observation model
# sd_1   - float  - standard deviation of likelihood 
# sd_2   - vector - standard deviation of prior
# N_e    - int    - number of samples to take
# logMar - bool   - whether to optimise the log posterior or marginal
fit.model_o = function(t,Y,W,sd_1,sd_2,N_e,logMar=F)
{
	## fit observation model
	ensemble = NULL
	for(k in 1:N_e)
	{
	    Omega_0 = rnorm(W*3,0,sd_2)
		if (logMar==F)
		{
	    	Omega_f   = argmax.logPost_o(t,Y,Omega_0,sd_1,sd_2) 
	    	logPost_0 = logPost_o(t,Y,Omega_0,sd_1,sd_2)
	    	logPost_f = logPost_o(t,Y,Omega_f,sd_1,sd_2)
		} else
		{
			Omega_f   = argmax.logMar_o(t,Y,Omega_0) 
	    	logPost_0 = logMar_o(t,Y,Omega_0)
	    	logPost_f = logMar_o(t,Y,Omega_f)
		}
	    ensemble = rbind(ensemble,c(logPost_f,Omega_f))
	    message(paste(k,"/",N_e,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
	}

	## terminate 
	return(ensemble)
}

## predict.model_o ##
## goal: compute predictions of the observation model over multiple time steps
## out:  returns list of predicted response and temporal derivative
# t        - vector - interpolated time steps
# t_       - vector - interpolated time steps
# ensemble - matrix - matrix containing multiple sampled parameter vectors
predict.model_o = function(t_,ensemble)
{
	## maximum a posteriori (MaP)
	idx            = which.max(ensemble[,1])
	MaP.Omega_o    = ensemble[idx,-1]
	MaP.Yhat_o     =     f_o.eval(t_,MaP.Omega_o)

	## compute response 
	Yhat_o         = t(apply(ensemble[,-1],1,function(x) f_o.eval(t_,x)))
	E.Yhat_o       = apply(Yhat_o,2,mean)
	q05.Yhat_o     = apply(Yhat_o,2,quantile,p=0.05)
	q95.Yhat_o     = apply(Yhat_o,2,quantile,p=0.95)
	
	## terminate
	results = list("MaP.Yhat_o"     = MaP.Yhat_o,
	 		 	   "E.Yhat_o"       = E.Yhat_o,
	 		 	   "q05.Yhat_o"     = q05.Yhat_o,
	 		 	   "q95.Yhat_o"     = q95.Yhat_o
		 		  )
	return(results)
}


## ddt.predict.model_o ##
## goal: compute predictions of the observation model over multiple time steps
## out:  returns list of predicted response and temporal derivative
# t        - vector - interpolated time steps
# t_       - vector - interpolated time steps
# ensemble - matrix - matrix containing multiple sampled parameter vectors
ddt.predict.model_o = function(t_,ensemble)
{
	## maximum a posteriori (MaP)
	idx            = which.max(ensemble[,1])
	MaP.Omega_o    = ensemble[idx,-1]
	MaP.ddt.Yhat_o = ddt.f_o.eval(t_,MaP.Omega_o)

	## compute derivative of response wtr to time
	ddt.Yhat_o     = t(apply(ensemble[,-1],1,function(x) ddt.f_o.eval(t_,x)))
	E.ddt.Yhat_o   = apply(ddt.Yhat_o,2,mean)
	q05.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.05)
	q95.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.95)

	## terminate
	results = list("MaP.ddt.Yhat_o" = MaP.ddt.Yhat_o,
	 		 	   "E.ddt.Yhat_o"   = E.ddt.Yhat_o,
	 		 	   "q05.ddt.Yhat_o" = q05.ddt.Yhat_o,
	 		       "q95.ddt.Yhat_o" = q95.ddt.Yhat_o
	 		 	  )
	return(results)
}

## plot.model_o ##
## goal: plot predictions of observation model
# t       - vector - time steps
# t_      - vector - interpolated time steps
# Y       - vector - observations 
# results - list   - predictions of observation model (i.e. interpolated response and temporal derivative)
# ...
# index   - vector - index to be added to the topright corner of each plot
#
## next:  - add main on top
##        - add legend
plot.model_o = function(t,t_,Y,results,col="red",xlab=c("Time","Time"),ylab=c("Y(t)","dY/dt"),index=NULL)
{
	## attach results
	attach(results,warn.conflicts=F)
	#
	## interpolated response
	plot(t_,rep(0,length(t_)),ylim=c(min(Y),max(Y))+0.2*c(-min(Y),max(Y)),type="l",lty=3,xlab=xlab[1],ylab=ylab[1])
	polygon(c(t_,rev(t_)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col,alpha=0.25),border=NA)
	lines(t_,E.Yhat_o,col=adjustcolor(col,0.75))
	points(t,Y)#,pch=16,col=adjustcolor("black",0.5))
	if(!is.null(index)) legend("topright",legend=index[1],bty="n",cex=1)
	#
	## temporal derivative
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*1.5,type="l",lty=3,xlab=xlab[2],ylab=ylab[2])
	polygon(c(t_,rev(t_)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col,alpha=0.25),border=NA)
	lines(t_,E.ddt.Yhat_o,col=adjustcolor(col,0.75))
	if(!is.null(index)) legend("topright",legend=index[2],bty="n",cex=1)
}

## DFT.NODEGM_o ##
## goal: compute interpolation/extrapolation and derivative of response variable of observation model using Discrete Fourier Transforms (DFT)
# t   - vector - time steps
# t_  - vector - interpolated/extrapolated time steps
# Y   - vector - values of response variables
# W   - int    - number of periodic elements in the series
# log - bool   - whether to log the response variable (e.g. if positive support only)
DFT.NODEGM_o = function(t,t_,Y,W,log=F)
{
	results = list()
	for(i in 1:ncol(Y))
	{
		DFT_         = fit.DFT(Y[,i],t,t_,W,log)
		results[[i]] = list("MaP.Yhat_o"     = DFT_[,1],
						    "MaP.ddt.Yhat_o" = DFT_[,2],	
						    "E.Yhat_o"       = DFT_[,1],
						    "E.ddt.Yhat_o"   = DFT_[,2],
						    "q05.Yhat_o"     = DFT_[,1],		
	                        "q05.ddt.Yhat_o" = DFT_[,2],
						    "q95.Yhat_o"     = DFT_[,1],		
	                        "q95.ddt.Yhat_o" = DFT_[,2])
	}
	return(results)
}

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## next:
##        - add switch between activation functions

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

## activation function of process model
f_sigma_p     = sigmo     # expo     # 
ddx.f_sigma_p = ddx.sigmo # ddu.expo # 

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
# Omega  - vector - vector of parameters of the network 
f_p.eval           = function(X,Omega) apply(t(X),2,function(x) f_p(x,Omega))
ddx.f_p.eval       = function(X,Omega) apply(t(X),2,function(x) ddx.f_p(x,Omega))
ddOmega.f_p.eval   = function(X,Omega) apply(t(X),2,function(x) ddOmega.f_p(x,Omega))
mean.ddx.f_p       = function(X,Omega) apply(ddx.f_p.eval(X,Omega),1,mean)

## Geber.f_p ##
## goal: compute contributions of each explanatory variable to the dynamics of the response at each time steps via Geber method 
# X      - matrix - matrix of variables across all time steps
# ddt.X  - matrix - matrix of temporal derivatives of variables
# Omega  - vector - vector of parameters of the network 
Geber.f_p = function(X,ddt.X,Omega) ddt.X * t(ddx.f_p.eval(X,Omega))

## prop.Geber.f_p ##
## goal: compute contributions as proportions of total contributions at each time step 
# X      - matrix - matrix of variables across all time steps
# ddt.X  - matrix - matrix of temporal derivatives of variables
# Omega  - vector - vector of parameters of the network 
prop.Geber.f_p = function(X,ddt.X,Omega) t(apply(Geber.f_p(X,ddt.X,Omega),1,function(x)(x^2)/sum(x^2)))

## prop.Geber.f_p ##
## goal: compute contributions as proportions of total contributions accross all time steps
# X      - matrix - matrix of variables across all time steps
# ddt.X  - matrix - matrix of temporal derivatives of variables
# Omega  - vector - vector of parameters of the network 
prop.tot.Geber.f_p = function(X,ddt.X,Omega)
{ 
	ssq = apply(Geber.f_p(X,ddt.X,Omega),2,function(x)sum(x^2))
	return(ssq/sum(ssq))
}

## r2_p ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
r2_p = function(X,Y,Omega)
{
  Yhat     = f_p.eval(X,Omega)
  res      = Y - Yhat
  r2       = 1 - sd(res^2)/sd(Y^2)
  return(r2)
}

## logLik_p ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
logLik_p = function(X,Y,Omega,sd_1)
{
  Yhat     = f_p.eval(X,Omega)
  res      = Y - Yhat
  logLik   = - sum((res^2)/(sd_1^2))
  return(logLik)
}

## logPrior_p ##
## goal: compute log prior density of the process model 
# Omega - vector - parameters
# sd_2  - float  - standard deviation of prior
logPrior_p = function(Omega,sd_2)
{
  logPrior = - sum((Omega^2)/(sd_2^2))
  return(logPrior)
}

## logPost_p ##
## goal: compute log posterior density of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
# sd_2  - float  - standard deviation of prior
logPost_p = function(X,Y,Omega,sd_1,sd_2)
{
  Yhat     = f_p.eval(X,Omega)
  res      = Y - Yhat
  logLik   = - sum((res^2)/(sd_1^2))
  logPrior = - sum((Omega^2)/(sd_2^2))
  logPost  = logLik + logPrior
  return(logPost)
}

## ddOmega.logPost_p ##
## goal: compute derivative of the log posterior density of the process model wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
# sd_2  - float  - standard deviation of prior
ddOmega.logPost_p = function(X,Y,Omega,sd_1,sd_2)
{
  Yhat             = f_p.eval(X,Omega)
  res              = Y - Yhat
  ddOmega.res      = - ddOmega.f_p.eval(X,Omega)
  ddOmega.logLik   = - 2 * ddOmega.res%*%res/(sd_1^2)
  ddOmega.logPrior = - 2 * Omega/(sd_2^2)
  ddOmega.logPost  = ddOmega.logLik + ddOmega.logPrior
  return(ddOmega.logPost)
}

## argmax.logPost_p ##
## goal: compute the maximum of the log posterior density of the process model
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
# sd_2  - float  - standard deviation of prior
argmax.logPost_p = function(X,Y,Omega,sd_1,sd_2)
{
	error_     = function(x) - logPost_p(X,Y,x,sd_1,sd_2)
	graderror_ = function(x) - ddOmega.logPost_p(X,Y,x,sd_1,sd_2)
	Omega      = optim(par    = Omega,
				       fn     = error_,
				       gr     = graderror_,
				       method = "BFGS"#,
				       # control=list("trace"=1,"REPORT"=1,"maxit"=100)
				       )$par
	return(Omega)
}

## logMar_p ##
## goal: compute log marginal posterior density of the process model
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
logMar_p = function(X,Y,Omega)
{
  Yhat     = f_p.eval(X,Omega)
  res      = Y - Yhat
  logLik   = - 0.5 * length(res)   * log(0.5 * sum(res^2) + 1)
  logPrior = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
  logMar   = logLik + logPrior
  return(logMar)
}

## ddOmega_p.logMar_p ##
## goal: compute derivative of the log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
ddOmega.logMar_p = function(X,Y,Omega)
{
	Yhat             = f_p.eval(X,Omega)
	res              = Y - Yhat 
	ddOmega.res      = - ddOmega.f_p.eval(X,Omega)
	ddOmega.logLik   = - 0.5 * length(res)   * 1/(0.5 * sum(res^2) + 1)   * 0.5 * ddOmega.res%*%res
	ddOmega.logPrior = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
	ddOmega.logMar   = ddOmega.logLik + ddOmega.logPrior
	return(ddOmega.logMar)
}

## argmax.logMar_p ##
## goal: compute the maximum of the log posterior density of the process model
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
argmax.logMar_p = function(X,Y,Omega)
{
	error_     = function(x) - logMar_p(X,Y,Omega)
	graderror_ = function(x) - ddOmega.logMar_p(X,Y,Omega)
	Omega      = optim(par    = Omega,
				       fn     = error_,
				       gr     = graderror_,
				       method = "BFGS"#,
				       # control=list("trace"=1,"REPORT"=1,"maxit"=100)
				       )$par
	return(Omega)
}

## fit.model_p ##
## goal: fit process model to dynamics of state variable (i.e. temporal derivative of interpolated time series)
## out:  returns list containing metadata regarding the scaling of variables and an ensemble of parameter vectors 
# X      - matrix - explanatory variables by columns (e.g. interpolated state variables)
# Y      - vector - response variable (e.g. temporal derivative of an interpolated state variable)
# W      - int    - number of hidden nodes in the process model
# sd_1   - float  - standard deviation of likelihood 
# sd_2   - float  - standard deviation of prior
# N_e    - int    - number of samples to take
fit.model_p = function(X,Y,W,sd_1,sd_2,N_e)
{
	## fit model
	ensemble= NULL 
	for(k in 1:N_e)    
	{   
	    Omega_0   = rnorm(W*(2+ncol(X))*2,0,sd_2)
	    Omega_f   = argmax.logPost_p(X,Y,Omega_0,sd_1,sd_2)
	    logPost_0 =        logPost_p(X,Y,Omega_0,sd_1,sd_2)
	    logPost_f =        logPost_p(X,Y,Omega_f,sd_1,sd_2)
	    ensemble  = rbind(ensemble,Omega_f)
	    message(paste(k,"/",N_e,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
	}
	## terminate
	return(ensemble)
}

## predict.model_p ##
## goal: predict response for values of X explanatory variables
# X        - matrix - explanatory variables (e.g. state variables)
# Y        - vector - response variable
# ensemble - matrix - list containing ensemble of parameter vectors
predict.model_p = function(X,ensemble)
{
	## compute predicted response 
	Yhat_p         = t(apply(ensemble,1,function(x) f_p.eval(X,x)))
	E.Yhat_p       = apply(Yhat_p,2,mean)
	q05.Yhat_p     = apply(Yhat_p,2,quantile,p=0.05)
	q95.Yhat_p     = apply(Yhat_p,2,quantile,p=0.95)

	## store
	predictions_p  = list("E.Yhat_p"   = E.Yhat_p, 
    		 	          "q05.Yhat_p" = q05.Yhat_p,    
    		 	          "q95.Yhat_p" = q95.Yhat_p 
    		 	         )
	## terminate
    return(predictions_p)
}

## ddx.predict.model_p ##
## goal: compute derivative of predicted response variable wtr to X explanatory variables
# X        - matrix - explanatory variables (e.g. state variables)
# ensemble - matrix - ensemble of parameter vectors
ddx.predict.model_p = function(X,ensemble)
{
	## compute derivative of predicted response wtr explanatory variables
	ddx.Yhat_p       = t(apply(ensemble,1,function(x) ddx.f_p.eval(X,x)))
	E.ddx.Yhat_p     = t(matrix(apply(ddx.Yhat_p,2,mean),ncol=nrow(X)))
	q05.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.05),ncol=nrow(X)))
	q95.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.95),ncol=nrow(X)))
		
	## store
	ddx.predictions_p = list("E.ddx.Yhat_p"    = E.ddx.Yhat_p,
    			            "q05.ddx.Yhat_p"   = q05.ddx.Yhat_p,
    			            "q95.ddx.Yhat_p"   = q95.ddx.Yhat_p
    			           )
	return(ddx.predictions_p)
}

## plot.model_p ##
## goal: visualise results of process model 
# t_     - vector - time steps
# Y      - vector - response variable that the process model has been fitted to
# ...
# index  - vector - index of graphs
# legend - vector - legend of different effects and contributions 
plot.model_p = function(t_,Y,predictions_p,col="red",xlab=c("Time","Time","Time"),ylab=c("Dynamics","Effects","Contributions"),index=NULL,legend=NULL)
{
	## initiate
	attach(predictions_p,warn.conflicts=F)
	N = ncol(E.ddx.Yhat_p)
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(Y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=ylab[1])
	points(t_,Y)#,pch=16,col=adjustcolor("black",0.5))
	polygon(c(t_,rev(t_)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_p,col=adjustcolor(col,alpha=0.75))
	if(!is.null(index)) legend("topright",legend=index[1],bty="n",cex=1)
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=ylab[2])
	for(j in 1:N) lines(t_,E.ddx.Yhat_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[2],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft",legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=ylab[3])
	for(j in 1:N) lines(t_, E.Geber_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[3],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft",legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
}

#
###
