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

## next:
##            - continue to clean and simplify code where possible (esp main functions)

#######################
## UTILITY FUNCTIONS ##
#######################

## goal: functions for visualisation

## .getmat ##
## goal: return matrix of function evaluated at grid points in variable limits 
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
.getMat = function(lims,func)
{
  ## compute 2D matrix
  res = 100
  X = apply(lims,1,function(X){seq(X[1],X[2],(X[2]-X[1])/res)}) # grid of pairwise combinations of x and y
  #
  ## fill output matrix - function value at each c(x,y)
  YMat = matrix(0,ncol=res+1,nrow=res+1)
  for(i in 1:length(X[,1]))
  { 
    for(j in 1:length(X[,2]))
    { 
      YMat[i,j] = func(c(X[i,1],X[j,2]))
    }
  }
  ## return matrix
  return(YMat)

}

## .heatmat ##
## goal: heatmap function for plotting 2D matrix
# func - function - function with two arguments
# lims - vector   - min and max of the two variables
# labs - vector   - labels of the x and y axis
# main - vector   - main title of the plot
.heatMat <- function(YMat,lims=rbind(c(0,1),c(0,1)),labs=c("",""),main=c(""),axes=T,maxAbsMinMax=NULL)
{
  ## plot matrix
  #
  ## relative min and max to be matched to min and max color levels
  if (is.null(maxAbsMinMax))  maxAbsMinMax <- max(abs(c(min(YMat),max(YMat))))
  #
  ## compute levels for color breaks
  levels = seq(-maxAbsMinMax,maxAbsMinMax,2*maxAbsMinMax/1000)
  colorLevels = rev(rainbow(1000,start=0,end=1,alpha=0.5))
  #
  ## heatmap
  image(YMat,breaks=levels,col=colorLevels,xaxt="n",yaxt="n",xlab=labs[1],ylab=labs[2],main=main)
  #
  ## add contour lines
  contour(YMat,add=T,col="black")
  #
  ## axes
  if(axes==T)
  { 
    for(i in 1:2){axis(side=i,at=c(0,.25,.5,.75,1),labels=round(c(0,.25,.5,.75,1)*diff(lims[i,])+lims[i,1],2))}
  }
}

## .plotMat ##
## goal: heatmap function for functions of two variables
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
# labs - vector   - labels of the x and y axis
# main - vector   - main title of the plot
.plotMat <- function(lims,func,labs=c("",""),main=c(""),axes=T,maxAbsMinMax=NULL)
{ 
  YMat = .getMat(lims,func)
  .heatMat(YMat,lims,labs,main,axes,maxAbsMinMax = maxAbsMinMax)
}

#
###

###################
## DFT FUNCTIONS ##
###################

## goal: functions to perform discrete fourier transforms

## DFT
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

## dDFT
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

## FFT
## goal: compute discrete Fourier transform of a signal f using the FFT algorithm
# f  - vector - vector of values of signal
# x  - vector - vector of dependent variable (e.g. time)
# x_ - vector - vector of dependent variable values at which to interpolate the signal
# K  - int    - number of elements in the series
FFT = function(f,x,x_,K)
{
    ## FFT
    dff = fft(f)

    ## upsample
    ndff = array(data = 0, dim = c(length(x_), 1L))
    ndff[1] = dff[1]
    ndff[2:(K+1)] = dff[2:(K+1)]
    ndff[length(ndff):(length(ndff) - K + 1)] = dff[length(f):(length(f) - K + 1)]

    ## frequency -> time domain
    indff   = fft(ndff/length(x), inverse = TRUE)
    return(indff)
}

## dFFT
## goal: compute the derivative of discrete Fourier transform of a signal f wtr to dependent variable x using the FFT algorithm
# f  - vector - vector of values of signal
# x  - vector - vector of dependent variable (e.g. time)
# x_ - vector - vector of dependent variable values at which to interpolate the signal
# K  - int    - number of elements in the series
dFFT = function(f,x,x_,K)
{
    ## FFT
    dff = fft(f)

    ## upsample
    ndff = array(data = 0, dim = c(length(x_), 1L))
    ndff[1] = dff[1]
    ndff[2:(K+1)] = dff[2:(K+1)]
    ndff[length(ndff):(length(ndff) - K + 1)] = dff[length(f):(length(f) - K + 1)]

    ## format kappa (# kappa = fftshift(kappa) # Re-order fft frequencies)
    kappa = (2*pi)*(1:length(x_))/max(x)
    m = length(kappa)
    p = ceiling(m/2)
    idx = c((p + 1):m, 1:p)
    kappa = kappa[idx]

    ## compute derivative
    dndff   = 1i*kappa*ndff

    ## frequency -> time domain
    idndff  = fft(dndff/length(x), inverse = TRUE)
    return(idndff)
}

## fit.DFT
## goal: 
fit.DFT = function(Y,x,x_,K,log=F)
{
	## standardise timesteps [min(t),max(t)] -> [1,N]
	t  = seq(1,length(x))
	nt = (x_ - min(x))/(max(x)-min(x))*(length(x)-1) + 1

	## standardise data
	if (log==T) Y_ = log(Y) else Y_ = Y
	mean_ = mean(Y_)
	sd_   = sd(Y_)
	Y_    = (Y_-mean_)/sd_
	
	## compute DFT
	Yhat_o_DFT     = DFT(Y_,t,nt,K)
	ddt.Yhat_o_DFT = dDFT(Y_,t,nt,K)
	
	## de-scale signal
	dt = diff(x[1:2])
	if (log==T)
	{	
		Yhat_o_DFT     = exp(mean_ + sd_ * Yhat_o_DFT)
		ddt.Yhat_o_DFT = 1/dt * sd_ * Yhat_o_DFT * ddt.Yhat_o_DFT
	} else
	{
		Yhat_o_DFT     = mean_ + sd_ * Yhat_o_DFT
		ddt.Yhat_o_DFT = 1/dt * sd_ * ddt.Yhat_o_DFT
	}

	## terminate
	return(cbind(Yhat_o_DFT,ddt.Yhat_o_DFT))
}


#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## goal: functions for the observation model

## f_o
## goal: interpolate state variable
## t     - float - time step in arbitrary units
## Omega - vect  - vector of parameters 
## output:
## float - value of the interpolation at time t
f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o
## goal: time derivative of the interpolated variable
## t     - float - time step in arbitrary units
## Omega - vect - vector of parameters 
## output:
## float - value of the derivative wtr to t of the interpolation at time t
ddt.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega_o.f_o
## goal: derivative of the interpolation wtr to each network parameter
## t     - float - time step in arbitrary units
## Omega - vect - vector of parameters 
## output:
## vect  - value of the derivative of the interpolation wtr to each parameter in the network
ddOmega_o.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	dfdOmega1 = sin(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
}

## evaluate functions across multiple time steps
## goal: evaluate functions across multiple time steps
## t     - float - time step in arbitrary units
## Omega - vect - vector of parameters 
## output:
## values of the corresponding functions evaluated across multiple time steps
f_o.eval           = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval       = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega_o.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega_o.f_o(x,Omega))

## logPost_o
## goal: evaluate the log posterior of the observation model given the observed data
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters
## output:
## float   - value of the log of the posterior of observation model
logPost_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	Yhat_o     = f_o.eval(t,Omega_o)
	res_o      = Y - Yhat_o
	logLik_o   = - sum((res_o^2)/(sd1_o^2))
	logPrior_o = - sum((Omega_o^2)/(sd2_o^2))
	logPost_o  = logLik_o + logPrior_o
	return(logPost_o)
}

## ddOmega_o.logPost_o
## goal: compute the derivate of the log posterior density wtr to each parameter
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters
## output:
## vector  - value of the derivative for each parameter
ddOmega_o.logPost_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	Yhat_o               = f_o.eval(t,Omega_o)
	res_o                = (Y - Yhat_o)
	ddOmega_o.res_o      = - ddOmega_o.f_o.eval(t,Omega_o)
	ddOmega_o.logLik_o   = - 2 * ddOmega_o.res_o%*%res_o/(sd1_o^2)
	ddOmega_o.logPrior_o = - 2 * Omega_o/(sd2_o^2)
	ddOmega_o.logPost_o  = ddOmega_o.logLik_o + ddOmega_o.logPrior_o
	return(ddOmega_o.logPost_o)
}

## argmax.logPost_o
## goal: maximise the log posterior density of the parameters
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - initial vector of parameters 
## output:
## Omega_o - vect  - vector of parameters that maximise (locally) the posterior
argmax.logPost_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	error_     = function(x) -logPost_o(t,Y,x,sd1_o,sd2_o)
	graderror_ = function(x) -ddOmega_o.logPost_o(t,Y,x,sd1_o,sd2_o)
	Omega_o  = optim(par    = Omega_o,
			 fn     = error_,
			 gr     = graderror_,
			 method = "BFGS"# ,
			 # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
			 )$par
	return(Omega_o)
}

## logMar_o
## goal: compute the log marginal posterior density of a parameter vector
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters
## output:
## float   - value of the log marginal posterior density
logMar_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	Yhat_o     = f_o.eval(t,Omega_o)
	res_o      = Y - Yhat_o
	logLik_o   = - 0.5 * length(Y) * log(0.5 * sum(res_o^2) + 1)
	logPrior_o = - 0.5 * length(Omega_o) * log(0.5 * sum(Omega_o^2) + 1)
	logMar_o   = logLik_o + logPrior_o
	return(logMar_o)
}

## ddOmega_o.logMar_o
## goal: compute the derivate of the log marginal posterior density wtr to each parameter
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters
## output:
## vector  - value of the derivative for eah parameter
ddOmega_o.logMar_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	Yhat_o               = f_o.eval(t,Omega_o)
	res_o                = (Y - Yhat_o)
	ddOmega_o.res_o      = -ddOmega_o.f_o.eval(t,Omega_o)
	ddOmega_o.logLik_o   = - 0.5 * length(Y) * 1/(0.5 * sum(res_o^2) + 1) * 0.5 * ddOmega_o.res_o%*%res_o
	ddOmega_o.logPrior_o = - 0.5 * length(Omega_o) * 1/(0.5 * sum(Omega_o^2) + 1) * Omega_o
	ddOmega_o.logMar_o   = ddOmega_o.logLik_o + ddOmega_o.logPrior_o
	return(ddOmega_o.logMar_o)
}

## argmax.logMar_o
## goal: maximise the log posterior density of the parameters
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - initial vector of parameters 
## output:
## Omega_o - vect  - vector of parameters that maximise (locally) the posterior
argmax.logMar_o = function(t,Y,Omega_o,sd1_o,sd2_o)
{
	error_     = function(x) -logMar_o(t,Y,x,sd1_o,sd2_o)
	graderror_ = function(x) -ddOmega_o.logMar_o(t,Y,x,sd1_o,sd2_o)
	Omega_o  = optim(par    = Omega_o,
			 fn     = error_,
			 gr     = graderror_,
			 method = "BFGS"# ,
			 # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
			 )$par
	return(Omega_o)
}
	
#
###

######################################
## MAIN FUNCTIONS OBSERVATION MODEL ##
######################################

## goal: functions to fit the observation and process model

## DFT.model_o
# DTF.model_o = function(Y,t,t_,K)

## fit.model_o
## goal: fit the observation model (i.e. interpolate the time series)
## t       - vector - vector of time steps
## Y       - vector - vector of state variable
## W_o     - int    - number of hidden nodes in the observation model
## sd1_o   - float  - standard deviation of gaussian likelihood of model prediction around observations 
## sd2_o   - float  - standard deviation of the prior distribution of model parameters 
## N_e     - int    - number of samples to take
## log     - bool   - whether to log the response variable (e.g. if its support is positive)
## logMar  - bool   - whether to optimise the log posterior or marginal
## output:
## Omega_o - vector - vector of parameters of the fitted observation model
fit.model_o = function(t,Y,W_o,sd1_o,sd2_o,N_e,log=F,logMar=F)
{
	## standardise time steps
	t = 1:length(t)

	## standardise/transform data
	if (log==T) Y = log(Y)
	mean_y = mean(Y)
	sd_y   = sd(Y)
	Y      = (Y-mean_y)/sd_y

	## fit observation model
	Omega_o_chain = NULL
	for(k in 1:N_e)
	{
	    Omega_o_0 = rnorm(W_o*3,0,sd2_o)
		if (logMar==F)
		{
	    	Omega_o_f   = argmax.logPost_o(t,Y,Omega_o_0,sd1_o,sd2_o) 
	    	logPost_o_0 = logPost_o(t,Y,Omega_o_0,sd1_o,sd2_o)
	    	logPost_o_f = logPost_o(t,Y,Omega_o_f,sd1_o,sd2_o)
		} else
		{
			Omega_o_f   = argmax.logMar_o(t,Y,Omega_o_0,sd1_o,sd2_o) 
	    	logPost_o_0 = logMar_o(t,Y,Omega_o_0,sd1_o,sd2_o)
	    	logPost_o_f = logMar_o(t,Y,Omega_o_f,sd1_o,sd2_o)
		}
	    Omega_o_chain   = rbind(Omega_o_chain,c(logPost_o_f,mean_y,sd_y,Omega_o_f))
	    message(paste(k,"/",N_e,"\t",
	        round(logPost_o_0,2),"\t","-->","\t",
	        round(logPost_o_f,2),sep=""))
	}

	## terminate 
	message("")
	message(paste("MaP: ",max(Omega_o_chain[,1],sep="")))
	message("")
	return(Omega_o_chain)
}

## predict.model_o
## goal: predict the observation model over multiple time steps
## t             - vector - vector of time steps
## Y             - vector - vector of state variable
## Omega_o_chain - matrix - matrix containing multiple sampled parameter vectors
## log           - bool   - whether to log the response variable (e.g. if its support is positive)
predict.model_o = function(t,Omega_o_chain,log=F)
{
	## standardise time steps
	t_real  = t
	dt_real = diff(t_real[1:2])
	t       = 1:length(t)
	
	## scaling
	mean_y = Omega_o_chain[1,2]
	sd_y   = Omega_o_chain[1,3]

	## maximum a posteriori (MaP)
	idx            = which.max(Omega_o_chain[,1])
	Omega_o_MaP    = Omega_o_chain[idx,-c(1:3)]
	MaP.Yhat_o     = f_o.eval(t,Omega_o_MaP)
	MaP.ddt.Yhat_o = ddt.f_o.eval(t,Omega_o_MaP)

	## compute response 
	Yhat_o         = t(apply(Omega_o_chain[,-c(1:3)],1,function(x) f_o.eval(t,x)))
	E.Yhat_o       = apply(Yhat_o,2,mean)
	q05.Yhat_o     = apply(Yhat_o,2,quantile,p=0.05)
	q95.Yhat_o     = apply(Yhat_o,2,quantile,p=0.95)
	
	## compute derivative of response wtr to time
	ddt.Yhat_o     = t(apply(Omega_o_chain[,-c(1:3)],1,function(x) ddt.f_o.eval(t,x)))
	E.ddt.Yhat_o   = apply(ddt.Yhat_o,2,mean)
	q05.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.05)
	q95.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.95)

	## de-scale data
	if(log==T)
	{
		MaP.Yhat_o     = exp(mean_y + sd_y * MaP.Yhat_o)
		E.Yhat_o       = exp(mean_y + sd_y * E.Yhat_o)
		q05.Yhat_o     = exp(mean_y + sd_y * q05.Yhat_o)
		q95.Yhat_o     = exp(mean_y + sd_y * q95.Yhat_o)
		MaP.ddt.Yhat_o = 1/dt_real * sd_y * MaP.Yhat_o * MaP.ddt.Yhat_o
		E.ddt.Yhat_o   = 1/dt_real * sd_y * MaP.Yhat_o * E.ddt.Yhat_o
		q05.ddt.Yhat_o = 1/dt_real * sd_y * MaP.Yhat_o * q05.ddt.Yhat_o
		q95.ddt.Yhat_o = 1/dt_real * sd_y * MaP.Yhat_o * q95.ddt.Yhat_o
	} else
	{
		MaP.Yhat_o     = mean_y + sd_y * MaP.Yhat_o
		E.Yhat_o       = mean_y + sd_y * E.Yhat_o
		q05.Yhat_o     = mean_y + sd_y * q05.Yhat_o
		q95.Yhat_o     = mean_y + sd_y * q95.Yhat_o
		MaP.ddt.Yhat_o = 1/dt_real * sd_y * MaP.ddt.Yhat_o
		E.ddt.Yhat_o   = 1/dt_real * sd_y * E.ddt.Yhat_o
		q05.ddt.Yhat_o = 1/dt_real * sd_y * q05.ddt.Yhat_o
		q95.ddt.Yhat_o = 1/dt_real * sd_y * q95.ddt.Yhat_o
	}

	## terminate
	pred_o        = list(MaP.Yhat_o,E.Yhat_o,q05.Yhat_o,q95.Yhat_o,MaP.ddt.Yhat_o,E.ddt.Yhat_o,q05.ddt.Yhat_o,q95.ddt.Yhat_o)
	names(pred_o) = c("MaP.Yhat_o","E.Yhat_o","q05.Yhat_o","q95.Yhat_o","MaP.ddt.Yhat_o","E.ddt.Yhat_o","q05.ddt.Yhat_o","q95.ddt.Yhat_o")
	return(pred_o)
}

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## activation functions
f_sigma     = function(x) 1/(1+exp(-x)) 
ddu.f_sigma = function(x) f_sigma(x) * (1 - f_sigma(x))
# f_sigma     = function(x) sin(2*pi*x) 
# ddu.f_sigma = function(x) 2*pi*cos(2*pi*x) 
# f_sigma     = function(x) x
# ddu.f_sigma = function(x) 1
# f_sigma     = function(x) x^2
# ddu.f_sigma = function(x) 2*x
# f_sigma     = function(x) x^3
# ddu.f_sigma = function(x) 3*x^2
# f_sigma     = function(x) (x>0)*x
# ddu.f_sigma = function(x) (x>0)
# f_sigma     = function(x) exp(-x^2)
# ddu.f_sigma = function(x) -2*x*exp(-x^2)

## f_p 
## goal: process model defined as a neural ordinary differential equation (NODE) 
##       to approximate dynamics of each state variable
## x     - vector - vector of input variables
## Omega - vector - vector of parameters of the NODE
## output:
## float - value of the NODE
f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2 + length(x))
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
}

## ddOmega_p.f_p
## goal: derivative of the process model wtr to parameters
## x      - vector - vector of input variables
## Omega  - vector - vector of parameters of the NODE
## output:
## vector - value of the derivative of the NODE wtr to each parameter
ddOmega_p.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	ddOmega_p1. = f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x)))
	ddOmega_p2. = Omega[,1]*1*ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x)))
	ddOmega_p3. = Omega[,1]*rep(1,nrow(Omega))%*%t(x)*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
	return(c(ddOmega_p1.,ddOmega_p2.,ddOmega_p3.))
}

## ddx.f_p
## goal: derivative of the process model wtr to state variables
## x      - vector - vector of input variables
## Omega  - vector - vector of parameters of the NODE
## output:
## vector - value of the derivative of the NODE wtr to each state variable
ddx.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	ddx. = Omega[,1]%*%(Omega[,-c(1:2)]*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x)))))
	return(ddx.)
}

## evaluate process functions across several timesteps
## X      - matrix - matrix of variables across all time steps
## Omega  - vector - vector of parameters of the network 
## output:
## vector - value of the function evaluated across all time steps
f_p.eval           = function(X,Omega) apply(t(X),2,function(x) f_p(x,Omega))
ddx.f_p.eval       = function(X,Omega) apply(t(X),2,function(x) ddx.f_p(x,Omega))
ddOmega_p.f_p.eval = function(X,Omega) apply(t(X),2,function(x) ddOmega_p.f_p(x,Omega))

## logPost_p
# goal: log posterior density of the process model (first level of inference)
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## sd1_p      - float  - standard deviation of process error
## sd2_p      - float  - standard deviation of prior
## output:
## float      - value of the posterior density givent the interpolated data and parameter vector
logPost_p = function(Yhat_o,ddt.Yhat_o,Omega_p,sd1_p,sd2_p)
{
  ddt.Yhat_p = f_p.eval(Yhat_o,Omega_p)
  res_p      = ddt.Yhat_o - ddt.Yhat_p
  logLik_p   = - sum((res_p^2)/(sd1_p^2))
  logPrior_p = - sum((Omega_p^2)/(sd2_p^2))
  logPost_p  = logLik_p + logPrior_p
  return(logPost_p)
}

## ddOmega_p.logPost_p
## goal: derivative of the log posterior density of the process model wtr to parameters
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## sd1_p      - float  - standard deviation of process error
## sd2_p      - float  - standard deviation of prior
## output:
## vector     - value of the derivative of the posterior density wtr to each parameter 
ddOmega_p.logPost_p = function(Yhat_o,ddt.Yhat_o,Omega_p,sd1_p,sd2_p)
{
  ddt.Yhat_p           = f_p.eval(Yhat_o,Omega_p)
  res_p                = ddt.Yhat_o - ddt.Yhat_p
  ddOmega_p.res_p      = - ddOmega_p.f_p.eval(Yhat_o,Omega_p)
  ddOmega_p.logLik_p   = - 2 * ddOmega_p.res_p%*%res_p/(sd1_p^2)
  ddOmega_p.logPrior_p = - 2 * Omega_p/(sd1_p^2)
  ddOmega_p.logPost_p  = ddOmega_p.logLik_p + ddOmega_p.logPrior_p
  return(ddOmega_p.logPost_p)
}

## argmax.logPost_p
## goal: find the maximum of the log posterior density of the process model
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## sd1_p      - float  - standard deviation of process error
## sd2_p      - float  - standard deviation of prior
## output:
## vector     - vector of parameters that maximise locally the posterior density
argmax.logPost_p = function(Yhat_o,ddt.Yhat_o,Omega_p,sd1_p,sd2_p)
{
	error_     = function(x) -logPost_p(Yhat_o,ddt.Yhat_o,x,sd1_p,sd2_p)
	graderror_ = function(x) -ddOmega_p.logPost_p(Yhat_o,ddt.Yhat_o,x,sd1_p,sd2_p)
	Omega_p  = optim(par=Omega_p,
				fn=error_,
				gr=graderror_,
				method="BFGS"#,
				# control=list("trace"=1,"REPORT"=1,"maxit"=100)
				)$par
	return(Omega_p)
}

## logMar_p
## goal: log marginal posterior density of the process model
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## output:
## float      - value of the log marginal posterior density given parameter vector and data
logMar_p = function(Yhat_o,ddt.Yhat_o,Omega_p,sd1_p,sd2_p)
{
  ddt.Yhat_p = f_p.eval(Yhat_o,Omega_p)
  res_p      = ddt.Yhat_o - ddt.Yhat_p
  logLik_p   = - 0.5 * length(ddt.Yhat_o) * log(0.5 * sum(res_p^2) + 1)
  logPrior_p = - 0.5 * length(Omega_p) * log(0.5 * sum(Omega_p^2) + 1)
  logMar_p   = logLik_p + logPrior_p
  return(logMar_p)
}

## ddOmega_p.logMar_p 
## goal: derivative of the log marginal posterior density wtr to parameters
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## output:
## vector     - vector - vector containing the derivative of the log marginal posterior density wtr to each parameter 
ddOmega_p.logMar_p = function(Yhat_o,ddt.Yhat_o,Omega_p,sd1_p,sd2_p)
{
	ddt.Yhat_p           = f_p.eval(Yhat_o,Omega_p)
	res_p                = ddt.Yhat_o - ddt.Yhat_p
	ddOmega_p.res_p      = - ddOmega_p.f_p.eval(Yhat_o,Omega_p)
	ddOmega_p.logLik_p   = - 0.5 * length(res_p) * 1/(0.5 * sum(res_p^2) + 1) * 0.5 * ddOmega_p.res_p%*%res_p
	ddOmega_p.logPrior_p = - 0.5 * length(Omega_p) * 1/(0.5 * sum(Omega_p^2) + 1) * Omega_p
	ddOmega_p.logMar_p   = ddOmega_p.logLik_p + ddOmega_p.logPrior_p
	return(ddOmega_p.logMar_p)
}

## argmax.logMar_p
## goal: find the maximum of the log posterior density of the process model
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
## sd1_p      - float  - standard deviation of process error
## sd2_p      - float  - standard deviation of prior
## output:
## vector     - vector of parameters that maximise locally the posterior density
argmax.logMar_p = function(Ybar_o,ddt.Ybar_o,Omega_p,sd1_p,sd2_p)
{
	error_     = function(x) -logMar_p(Ybar_o,ddt.Ybar_o,x,sd1_p,sd2_p)
	graderror_ = function(x) -ddOmega_p.logMar_p(Ybar_o,ddt.Ybar_o,x,sd1_p,sd2_p)
	Omega_p  = optim(par=Omega_p,
				fn=error_,
				gr=graderror_,
				method="BFGS"#,
				# control=list("trace"=1,"REPORT"=1,"maxit"=100)
				)$par
	return(Omega_p)
}

#
###

##################################
## MAIN FUNCTIONS PROCESS MODEL ##
##################################

## goal:

## fit.model_p
## goal: fit process model to dynamics of state variable (i.e. temporal derivative of interpolated time series)
## X      - matrix - table containing explanatory variables by columns (e.g. interpolated state variables)
## Y      - vector - table containing response variable (e.g. temporal derivative of an interpolated state variable)
## W_p    - int    - number of hidden nodes in the process model
## sd1_p  - float  - standard deviation of gaussian likelihood of process model prediction around observations 
## sd2_p  - float  - standard deviation of the prior distribution of model parameters 
## N_e    - int    - number of samples to take
## log    - bool   - whether to log the response variable (e.g. if its support is positive)
## logMar - bool   - whether to optimise the log posterior or marginal
## output:
## Omega_p    - vector - vector of parameters of the fitted process model
fit.model_p = function(X,Y,W_p,sd1_p,sd2_p,N_e,log=F,logMar=F)
{
	## standardise predictive variables
	X_     = X
	if (log == T) X_ = log(X_)
	mean_x = apply(X_,2,mean)
	sd_x   = apply(X_,2,sd)
	X_     = t((t(X_)-mean_x)/sd_x)
	
	## standardise response variable
	Y_     = Y
	mean_y = mean(Y_)
	sd_y   = sd(Y_)
	Y_     = (Y_-mean_y)/sd_y
	
	## fit model
	Omega_p_chain = NULL 
	for(k in 1:N_e)    
	{   
	    Omega_p_0 = rnorm(W_p*(2+ncol(X_)),0,sd2_p)
		if(logMar==F)
		{
	    	Omega_p_f      = argmax.logPost_p(X_,Y_,Omega_p_0,sd1_p,sd2_p)
	    	logPost_p_0    = logPost_p(X_,Y_,Omega_p_0,sd1_p,sd2_p)
	    	logPost_p_f    = logPost_p(X_,Y_,Omega_p_f,sd1_p,sd2_p)
		} else
		{
			Omega_p_f      = argmax.logMar_p(X_,Y_,Omega_p_0,sd1_p,sd2_p)
	    	logPost_p_0    = logMar_p(X_,Y_,Omega_p_0,sd1_p,sd2_p)
	    	logPost_p_f    = logMar_p(X_,Y_,Omega_p_f,sd1_p,sd2_p)
		}
	    Omega_p_chain  = rbind(Omega_p_chain,c(logPost_p_f,mean_y,sd_y,Omega_p_f))
	    message(paste(k,"/",N_e,"\t",
	            round(logPost_p_0,2),"\t","-->","\t",
	            round(logPost_p_f,2),sep=""))
	}
	
	## terminate
	message("") 
	message(max(Omega_p_chain[,1]))
	message("") 
	return(Omega_p_chain)
}

## predict.model_p
## goal: evaluate process model for values of X explanatory variables
## X             - matrix - matrix of explanatory variables (e.g. state variables)
## Y             - vector - vector containing values of response variable
## Omega_p_chain - matrix - matrix containing parameter vector sampled by row
## output:
predict.model_p = function(X,Omega_p_chain)
{
	## standardise predictive variables
	mean_x = apply(X,2,mean)
	sd_x   = apply(X,2,sd)
	X_     = t((t(X)-mean_x)/sd_x)
	
	## scaling
	mean_y = Omega_p_chain[1,2]
	sd_y   = Omega_p_chain[1,3]

	## maximum a posteriori
	idx         = which.max(Omega_p_chain[,1])
	Omega_p_MaP = Omega_p_chain[idx,-c(1:3)]
	
	## evaluate best model
	MaP.Yhat_p     = f_p.eval(X_,Omega_p_MaP)
	MaP.ddx.Yhat_p = t(ddx.f_p.eval(X_,Omega_p_MaP))
	
	## compute predicted response 
	Yhat_p         = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) f_p.eval(X_,x)))
	E.Yhat_p       = apply(Yhat_p,2,mean)
	q05.Yhat_p     = apply(Yhat_p,2,quantile,p=0.05)
	q95.Yhat_p     = apply(Yhat_p,2,quantile,p=0.95)

	## compute derivative of predicted response wtr explanatory variables
	ddx.Yhat_p     = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) ddx.f_p.eval(X_,x)))
	E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,mean),ncol=length(E.Yhat_p)))
	q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.05),ncol=length(E.Yhat_p)))
	q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.95),ncol=length(E.Yhat_p)))
	
	## de-scale predictions
	MaP.Yhat_p     = mean_y + sd_y * MaP.Yhat_p
	E.Yhat_p       = mean_y + sd_y * E.Yhat_p
	q05.Yhat_p     = mean_y + sd_y * q05.Yhat_p
	q95.Yhat_p     = mean_y + sd_y * q95.Yhat_p
	MaP.ddx.Yhat_p = sd_y * t(1/sd_x * t(MaP.ddx.Yhat_p))
	E.ddx.Yhat_p   = sd_y * t(1/sd_x * t(E.ddx.Yhat_p))
	q05.ddx.Yhat_p = sd_y * t(1/sd_x * t(q05.ddx.Yhat_p))
	q95.ddx.Yhat_p = sd_y * t(1/sd_x * t(q95.ddx.Yhat_p))
	
	## store
	pred_p        = list(E.Yhat_p,q05.Yhat_p,q95.Yhat_p,E.ddx.Yhat_p,q05.ddx.Yhat_p,q95.ddx.Yhat_p)
	names(pred_p) = c("E.Yhat_p","q05.Yhat_p","q95.Yhat_p","E.ddx.Yhat_p","q05.ddx.Yhat_p","q95.ddx.Yhat_p")
	return(pred_p)
}

## Geber.model_p
## goal: 
## input:
## output:
Geber.model_p = function(X,ddt.X=NULL,Omega_p_chain)
{
	## standardise predictive variables
	mean_x = apply(X,2,mean)
	sd_x   = apply(X,2,sd)
	X_     = t((t(X)-mean_x)/sd_x)

	## scaling
	mean_y = Omega_p_chain[1,2]
	sd_y   = Omega_p_chain[1,3]
	
	## maximum a posteriori
	idx         = which.max(Omega_p_chain[,1])
	Omega_p_MaP = Omega_p_chain[idx,-c(1:3)]
	
	## evaluate best model
	MaP.Yhat_p     = f_p.eval(X_,Omega_p_MaP)
	MaP.ddx.Yhat_p = t(ddx.f_p.eval(X_,Omega_p_MaP))
	
	## compute predicted response 
	Yhat_p         = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) f_p.eval(X_,x)))
	E.Yhat_p       = apply(Yhat_p,2,mean)
	q05.Yhat_p     = apply(Yhat_p,2,quantile,p=0.05)
	q95.Yhat_p     = apply(Yhat_p,2,quantile,p=0.95)
	
	## compute derivative of predicted response wtr explanatory variables
	ddx.Yhat_p     = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) ddx.f_p.eval(X_,x)))
	E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,mean),ncol=length(E.Yhat_p)))
	q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.05),ncol=length(E.Yhat_p)))
	q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.95),ncol=length(E.Yhat_p)))

	## compute geber
	if (is.null(ddt.X)) ddt.X = Yhat_p 
	Geber_p        = t(as.vector(t(ddt.X))*t(ddx.Yhat_p))
	E.Geber_p      = t(matrix(apply(Geber_p,2,mean),ncol=length(E.Yhat_p)))
	q05.Geber_p    = t(matrix(apply(Geber_p,2,quantile,p=0.05),ncol=length(E.Yhat_p)))
	q95.Geber_p    = t(matrix(apply(Geber_p,2,quantile,p=0.95),ncol=length(E.Yhat_p)))
	
	## de-scale predictions
	E.Geber_p      = sd_y * t(1/sd_x * t(E.Geber_p))
	q05.Geber_p    = sd_y * t(1/sd_x * t(q05.Geber_p))
	q95.Geber_p    = sd_y * t(1/sd_x * t(q95.Geber_p))
	
	## store
	pred_p         = list(E.Geber_p,q05.Geber_p,q95.Geber_p)
	names(pred_p)  = c("E.Geber_p","q05.Geber_p","q95.Geber_p")
	return(pred_p)
}

#
###

####################
## MAIN FUNCTIONS ##
####################

## goal: main functions to fit NODEs with gradient matching (GM)

## DFT.NODEGM_o
DFT.NODEGM_o = function(t,t_,Y,W_o,log=F)
{
	pred_o = list()
	for(i in 1:ncol(Y))
	{
		DFT_        = fit.DFT(Y[,i],t,t_,W_o,log)
		pred_o[[i]] = list("MaP.Yhat_o"     = DFT_[,1],
						   "MaP.ddt.Yhat_o" = DFT_[,2],	
						   "E.Yhat_o"       = DFT_[,1],
						   "E.ddt.Yhat_o"   = DFT_[,2],
						   "q05.Yhat_o"     = DFT_[,1],		
	                       "q05.ddt.Yhat_o" = DFT_[,2],
						   "q95.Yhat_o"     = DFT_[,1],		
	                       "q95.ddt.Yhat_o" = DFT_[,2])
	}
	return(pred_o)
}

## fit.NODEGM_o
fit.NODEGM_o = function(t,t_,Y,W_o,sd1_o,sd2_o,N_e,log=F,logMar=F)
{
	pred_o = list()
	for(i in 1:ncol(Y))
	{
		Omega_o_chain = fit.model_o(t,Y[,i],W_o,sd1_o,sd2_o,N_e,log,logMar)
		pred_o[[i]]   = predict.model_o(t_,Omega_o_chain,log)
	}
	return(pred_o)
}

## fit.NODEGM_p
fit.NODEGM_p = function(X,ddt.X,Y,W_p,sd1_p,sd2_p,N_e,log=F,logMar=F)
{
	pred_p = list()
	for(i in 1:ncol(Y))
	{
		Omega_p_chain = fit.model_p(X,Y[,i],W_p,sd1_p,sd2_p,N_e,log,logMar)
		pred_p[[i]]   = c(predict.model_p(X,Omega_p_chain),Geber.model_p(X,ddt.X,Omega_p_chain))
	}
	return(pred_p)
}

## plot.NODEGM_o
## goal:
plot.NODEGM_o = function(t,t_,Y,pred_o,col="red")
{
	attach(pred_o,warn.conflicts=F)
	#
	## interpolation
	plot(t_,rep(0,length(t_)),ylim=c(0,1)*max(abs(Y))*1.5,type="l",lty=3,ylab="State")
	points(t,Y,pch=16)
	polygon(c(t_,rev(t_)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_o,col=col)
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*1.5,type="l",lty=3,ylab="Dynamics")
	polygon(c(t_,rev(t_)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.ddt.Yhat_o,col=col)
}


## plot.NODEGM_p
## goal: 
plot.NODEGM_p = function(t,t_,Y,pred_p,col="red")
{
	attach(pred_p,warn.conflicts=F)
	N = ncol(E.ddx.Yhat_p)
	# par(mfrow=c(3,1),mar=c(2,5,1,1),cex.lab=1.5)
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(Y))*1.5,type="l",lty=3,ylab="Dynamics")
	points(t_,Y,pch=16)
	polygon(c(t_,rev(t_)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_p,col=col)
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,ylab="Effect")
	for(j in 1:N) lines(t_, E.ddx.Yhat_p[,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,ylab="Contribution")
	for(j in 1:N) lines(t_, E.Geber_p[,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	#
	# par(mfrow=c(1,1))
}

# ## plot.NODEGM
# ## goal: visualise the results of the NODEGM fit
# ## input:
# ## TS      - matrix - matrix containing the original time series data 
# ## results - list   - output list of the summary function containing the quantities to visualise (e.g. mean, q05, and q95 of effects)
# plot.NODEGM_o = function(t,t_,Y,results)
# {
# 	attach(results)
# 	N = ncol(Y)
# 	# t_ = t_/alpha_o*(max(t)-min(t)) + min(t)
# 	for(i in 1:N)
# 	{
# 		par(mfrow=c(5,1),mar=c(2,5,1,1),cex.lab=1.5)
# 		#
# 		## interpolation
# 		plot(t_,rep(0,length(t_)),ylim=c(-1,1)*3,type="l",lty=2,ylab="State")
# 		points(t,Y[,i])
# 		lines(t_,E.Ybar_o[,i],col=rainbow(N)[i])
# 		polygon(c(t_,rev(t_)),c(q05.Ybar_o[,i],rev(q95.Ybar_o[,i])),col=rainbow(N,alpha=0.2)[i],border=NA)
# 		#
# 		## dynamics
# 		plot(t_,rep(0,length(t_)),ylim=c(-1,1)*3,type="l",lty=2,ylab="Dynamics")
# 		lines(t_,1/E.Ybar_o[,i]*E.ddt.Ybar_o[,i],col="black")
# 		lines(t_,E.ddt.Ybar_p[,i],col=rainbow(N)[i])
# 		polygon(c(t_,rev(t_)),c(q05.ddt.Ybar_p[,i],rev(q95.ddt.Ybar_p[,i])),col=rainbow(N,alpha=0.2)[i],border=NA)
# 		#
# 		## effects
# 		plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2,ylab="Effect")
# 		for(j in 1:N) lines(t_, E.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
# 		for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],rev(q95.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
# 		#
# 		## Geber 
# 		plot(t_,rep(0,length(t_)),ylim=c(-1,1)*0.5,type="l",lty=2,ylab="Contribution")
# 		for(j in 1:N) lines(t_, E.Geber_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
# 		for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,(i-1)*N+(1:N)][,j],rev(q95.Geber_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
# 		#
# 		## relative Geber 
# 		plot(t_,rep(0,length(t_)),ylim=c(0,1),type="l",lty=2,ylab="Relative cont.")
# 		for(j in 1:N) lines(t_, E.prop.Geber_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
# 		for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.prop.Geber_p[,(i-1)*N+(1:N)][,j],rev(q95.prop.Geber_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
# 		#
# 		par(mfrow=c(1,1))
# 	}
# }

#
###
