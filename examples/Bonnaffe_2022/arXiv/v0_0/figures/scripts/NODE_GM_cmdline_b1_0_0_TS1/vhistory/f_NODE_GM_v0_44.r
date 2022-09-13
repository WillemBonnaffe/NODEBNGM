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

#######################
## UTILITY FUNCTIONS ##
#######################

## goal: utility functions

## util.getmat ##
## goal: return matrix of function values obtained at grid points in variable limits 
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
# res  - int      - size of the matrix
util.getMat = function(lims,func,res=100)
{
 	## compute 2D matrix
	X = apply(lims,1,function(x)seq(x[1],x[2],(x[2]-x[1])/res)) # grid of pairwise combinations of x and y

	## fill output matrix - function value at each c(x,y)
	YMat = matrix(0,ncol=res+1,nrow=res+1)
	for(i in 1:length(X[,1]))
	{ 
		for(j in 1:length(X[,2]))
		{ 
			YMat[i,j] = func(c(X[i,1],X[j,2]))
		}
	}

  ## terminate
  return(YMat)
}

## util.heatmat ##
## goal: plot a heatmap of a 2D matrix 
# func - function - function with two arguments
# lims - vector   - min and max of the two variables
# labs - vector   - labels of the x- and y-axis
# main - vector   - main title of the plot
util.heatMat = function(YMat,lims=rbind(0:1,0:1),labs=c("",""),main=c(""),axes=T,maxAbsMinMax=NULL)
{
  ## relative min and max to be matched to min and max color levels
  if (is.null(maxAbsMinMax))  maxAbsMinMax = max(abs(c(min(YMat),max(YMat))))
  
  ## compute levels for color breaks
  levels      = seq(-maxAbsMinMax,maxAbsMinMax,2*maxAbsMinMax/1000)
  colorLevels = rev(rainbow(1000,start=0,end=1,alpha=0.5))
  
  ## heatmap
  image(YMat,breaks=levels,col=colorLevels,xaxt="n",yaxt="n",xlab=labs[1],ylab=labs[2],main=main)
  
  ## add contour lines
  contour(YMat,add=T)
  
  ## axes
  if(axes==T) for(i in 1:2) axis(side=i,at=c(0,.25,.5,.75,1),labels=round(c(0,.25,.5,.75,1)*diff(lims[i,])+lims[i,1],2))
}

## util.plotMat ##
## goal: plot a heatmap of a function taking two variables as input
# lims         - matrix   - limits of the two variables listed by rows
# func         - function - function with two arguments
# labs         - vector   - labels of the x- and y-axis
# main         - vector   - main title of the plot
# axes         - bool     - whether to add axes to the plot
# maxAbsMinMax - float    - maximum absolute value to colour
util.plotMat <- function(lims,func,labs=c("",""),main=c(""),axes=T,maxAbsMinMax=NULL)
{ 
	YMat = util.getMat(lims,func)
	util.heatMat(YMat,lims,labs,main,axes,maxAbsMinMax = maxAbsMinMax)
}

## util.diamondplot ##
## goal: plot distributions of variables sorted by columns in a table
# df      - matrix - matrix containing variables to be plotted sorted by columns
# y_lim   - vector - vector containing min and max value of the x-axis
# colVect - vector - vector of colours
# y_lab   - string - label of the y-axis
# main    - string - title of the plot
util.plot.diamond = function(df,y_lim=NULL,colVect=NULL,y_lab="",main="")
{
	## initiate graphical parameters
	alpha = 0.2 # scaling factor for the plot
	n_col = ncol(df)
	n_row = nrow(df)
	x     = 1:n_col
	x_lim = c(min(x)-alpha,max(x)+alpha)
	tmp   = as.vector(apply(df,2,function(x)density(x)$x))
	if(is.null(y_lim))     y_lim = c(min(tmp),max(tmp));
	if(is.null(colVect)) colVect = rainbow(n_col,alpha=0.5) 
	y     = round(seq(y_lim[1],y_lim[2],diff(y_lim)/5),2)
  
	## plot
	plot(rep(1,n_row),df[,1],xlim=x_lim,ylim=y_lim,cex=0,xlab="",ylab=y_lab,yaxt="n",main=main)
	axis(side=1,at=x,labels=x,tick=T,line=0)
	for(j in 1:n_col)
	{
		## coordinates
		x_ = rep(j,n_row)
		y_ = df[,j]
		x_density = density(df[,j])$x
		y_density = density(df[,j])$y
		y_density = y_density/max(y_density)*alpha # normalise
		
		## draw distributions
		polygon(x=c(j-y_density,rev(j+y_density)),y=c(x_density,rev(x_density)),col=colVect[j],border=NA)
	}
	lines(0:(n_col+1),rep(0,n_col+2),lty=2)
}

## util.plot.Hdiamond ##
## goal: plot horizontally distributions of the variables in a table 
# df      - matrix - matrix containing variables to be plotted sorted by columns
# y_lim   - vector - vector containing min and max value of the x-axis
# colVect - vector - vector of colours
# y_lab   - string - label of the y-axis
# main    - string - title of the plot
util.plot.Hdiamond = function(df,y_lim=NULL,colVect=NULL,y_lab="",main="",labels=NULL)
{
	## reverse columns
	df = df[,ncol(df):1]
	
	## initiate graphical parameters
	alpha = 0.2 # scaling factor for the plot
	n_col = ncol(df)
	n_row = nrow(df)
	x     = 1:n_col
	x_lim = c(min(x)-alpha,max(x)+alpha)
	tmp   = as.vector(apply(df,2,function(x)density(x)$x))
	if(is.null(y_lim)) y_lim = c(min(tmp),max(tmp));
	if(is.null(colVect)) colVect = rainbow(n_col,alpha=0.5) 
	if(is.null(labels)) labels=x
	colVect = rev(colVect)
	labels  = rev(labels)
  
	## plot
	plot(df[,1],rep(1,n_row),xlim=y_lim,ylim=x_lim,cex=0,xlab="",ylab=y_lab,main=main,yaxt="n")
	axis(side=2,at=x,labels=labels,tick=T,line=0,las=1)
	for(j in 1:n_col)
	{
		## coordinates
		x_ = rep(j,n_row)
		y_ = df[,j]
		x_density = density(df[,j])$x
		y_density = density(df[,j])$y
		y_density = y_density/max(y_density)*alpha # normalise
		
		## draw distributions
		polygon(y=c(j-y_density,rev(j+y_density)),x=c(x_density,rev(x_density)),col=colVect[j],border=NA)
	}
	lines(rep(0,n_col+2),0:(n_col+1),lty=2)
}

#
###

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

## FFT ##
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
    ndff          = array(data = 0, dim = c(length(x_), 1L))
    ndff[1]       = dff[1]
    ndff[2:(K+1)] = dff[2:(K+1)]
    ndff[length(ndff):(length(ndff) - K + 1)] = dff[length(f):(length(f) - K + 1)]

    ## frequency -> time domain
    indff   = fft(ndff/length(x), inverse = TRUE)
    return(indff)
}

## dFFT ##
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
## goal: compute parameter vecctor that maximises the log marginal density
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
# log    - bool   - whether to log the response variable (e.g. if its support is positive)
# logMar - bool   - whether to optimise the log posterior or marginal
# 
## next:
#        - instead of embedding scaling in ensemble return as metadata of the ensemble (along with time steps)
fit.model_o = function(t,Y,W,sd_1,sd_2,N_e,log=F,logMar=F)
{
	## standardise time steps
	t_0 = min(t) 
	t_f = max(t)
	dt  = diff(t[1:2])
	t   = 1:length(t)

	## standardise data
	if (log==T) Y = log(Y)
	mean_y = mean(Y)
	sd_y   = sd(Y)
	Y      = (Y-mean_y)/sd_y

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
	            round(logPost_0,2),"\t","-->","\t",
	            round(logPost_f,2),sep=""))
	}

	## terminate 
	message("")
	message(paste("MaP: ",max(ensemble[,1],sep="")))
	message("")
	results = list("metaData" = list("t_0"=t_0,"t_f"=t_f,"dt"=dt,"mean_y"=mean_y,"sd_y"=sd_y,"log"=log),
				   "ensemble" = ensemble
				  )
	return(results)
}

## predict.model_o ##
## goal: compute predictions of the observation model over multiple time steps
## out:  returns list of predicted response and temporal derivative
# t        - vector - interpolated time steps
# t_       - vector - interpolated time steps
# ensemble - matrix - matrix containing multiple sampled parameter vectors
# log      - bool   - whether to log the response variable (e.g. if its support is positive)
predict.model_o = function(t_,ensembleList)
{
	## unpack
	metaData = ensembleList$metaData
	ensemble = ensembleList$ensemble

	## standardisation parameters
	t_0    = metaData[["t_0"]]
	t_f    = metaData[["t_f"]]
	dt     = metaData[["dt"]]
	mean_y = metaData[["mean_y"]]
	sd_y   = metaData[["sd_y"]]
	log    = metaData[["log"]]

	## standardise time steps [min(t),max(t)] -> [1,N_t]
    t_ = (t_-t_0)/(t_f-t_0)*((t_f-t_0)/dt) + 1	

	## maximum a posteriori (MaP)
	idx            = which.max(ensemble[,1])
	MaP.Omega_o    = ensemble[idx,-1]
	MaP.Yhat_o     =     f_o.eval(t_,MaP.Omega_o)
	MaP.ddt.Yhat_o = ddt.f_o.eval(t_,MaP.Omega_o)

	## compute response 
	Yhat_o         = t(apply(ensemble[,-1],1,function(x) f_o.eval(t_,x)))
	E.Yhat_o       = apply(Yhat_o,2,mean)
	q05.Yhat_o     = apply(Yhat_o,2,quantile,p=0.05)
	q95.Yhat_o     = apply(Yhat_o,2,quantile,p=0.95)
	
	## compute derivative of response wtr to time
	ddt.Yhat_o     = t(apply(ensemble[,-1],1,function(x) ddt.f_o.eval(t_,x)))
	E.ddt.Yhat_o   = apply(ddt.Yhat_o,2,mean)
	q05.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.05)
	q95.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.95)

	## de-standardise
	if(log==T)
	{
		MaP.Yhat_o     = exp(mean_y + sd_y * MaP.Yhat_o)
		E.Yhat_o       = exp(mean_y + sd_y * E.Yhat_o)
		q05.Yhat_o     = exp(mean_y + sd_y * q05.Yhat_o)
		q95.Yhat_o     = exp(mean_y + sd_y * q95.Yhat_o)
		MaP.ddt.Yhat_o = 1/dt * sd_y * MaP.Yhat_o * MaP.ddt.Yhat_o
		E.ddt.Yhat_o   = 1/dt * sd_y * MaP.Yhat_o * E.ddt.Yhat_o
		q05.ddt.Yhat_o = 1/dt * sd_y * MaP.Yhat_o * q05.ddt.Yhat_o
		q95.ddt.Yhat_o = 1/dt * sd_y * MaP.Yhat_o * q95.ddt.Yhat_o
	} else
	{
		MaP.Yhat_o     = mean_y + sd_y * MaP.Yhat_o
		E.Yhat_o       = mean_y + sd_y * E.Yhat_o
		q05.Yhat_o     = mean_y + sd_y * q05.Yhat_o
		q95.Yhat_o     = mean_y + sd_y * q95.Yhat_o
		MaP.ddt.Yhat_o =   1/dt * sd_y * MaP.ddt.Yhat_o
		E.ddt.Yhat_o   =   1/dt * sd_y * E.ddt.Yhat_o
		q05.ddt.Yhat_o =   1/dt * sd_y * q05.ddt.Yhat_o
		q95.ddt.Yhat_o =   1/dt * sd_y * q95.ddt.Yhat_o
	}

	## terminate
	results = list("MaP.Yhat_o"     = MaP.Yhat_o,
	 		 	   "E.Yhat_o"       = E.Yhat_o,
	 		 	   "q05.Yhat_o"     = q05.Yhat_o,
	 		 	   "q95.Yhat_o"     = q95.Yhat_o,
	 		 	   "MaP.ddt.Yhat_o" = MaP.ddt.Yhat_o,
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
	plot(t_,rep(0,length(t_)),ylim=c(0,1)*max(abs(Y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=ylab[1])
	polygon(c(t_,rev(t_)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_o,col=col)
	points(t,Y,pch=16,col=adjustcolor("black",0.5))
	if(!is.null(index)) legend("topright",legend=index[1],bty="n",cex=1)
	#
	## temporal derivative
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*1.5,type="l",lty=3,xlab=xlab[2],ylab=ylab[2])
	polygon(c(t_,rev(t_)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.ddt.Yhat_o,col=col)
	if(!is.null(index)) legend("topright",legend=index[2],bty="n",cex=1)
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

## f_p ##
## goal: process model defined as single layer perceptron with input x and parameter vector Omega 
# x     - vector - input variables
# Omega - vector - parameters
f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2 + length(x))
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
}

## ddOmega.f_p ##
## goal: compute derivative vector of the process model wtr to each parameter
# x      - vector - input variables
# Omega  - vector - parameters 
ddOmega.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	x       = t(x)
	Omega_1 = Omega[,1]
	Omega_2 = Omega[,2]
	Omega_3 = Omega[,-c(1:2)]
	ddOmega_1 = f_sigma(Omega_2 + Omega_3%*%t(x))
	ddOmega_2 = Omega_1 * ddu.f_sigma(Omega_2 + Omega_3%*%t(x))
	ddOmega_3 = Omega_1%*%x * as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x)))
	# ddOmega_3 = Omega[,1]*rep(1,nrow(Omega))%*%t(x)*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
	return(c(ddOmega_1,ddOmega_2,ddOmega_3))
}

## ddx.f_p ##
## goal: compute derivative vector of the process model wtr to each state variable
# x      - vector - input variables
# Omega  - vector - parameters
ddx.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	ddx = Omega[,1]%*%(Omega[,-c(1:2)]*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x)))))
	return(ddx)
}

## functions to compute across several timesteps
# X      - matrix - matrix of variables across all time steps
# Omega  - vector - vector of parameters of the network 
f_p.eval           = function(X,Omega)       apply(t(X),2,function(x) f_p(x,Omega))
ddx.f_p.eval       = function(X,Omega)       apply(t(X),2,function(x) ddx.f_p(x,Omega))
ddOmega.f_p.eval   = function(X,Omega)       apply(t(X),2,function(x) ddOmega.f_p(x,Omega))
mean.ddx.f_p       = function(X,Omega)       apply(ddx.f_p.eval(X,Omega),1,mean)

## Geber functions
Geber.f_p          = function(X,ddt.X,Omega) ddt.X * t(ddx.f_p.eval(X,Omega))
prop.Geber.f_p     = function(X,ddt.X,Omega) t(apply(Geber.f_p(X,ddt.X,Omega),1,function(x)(x^2)/sum(x^2)))
prop.tot.Geber.f_p = function(X,ddt.X,Omega)
{ 
	ssq = apply(Geber.f_p(X,ddt.X,Omega),2,function(x)sum(x^2))
	return(ssq/sum(ssq))
}

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
summary.model_p = function(X,ddt.X=NULL,Omega_p_chain)
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
	E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,mean),ncol=nrow(X_)))
	q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.95),ncol=nrow(X_)))

	## compute Geber
	if (is.null(ddt.X)) ddt.X = Yhat_p 
	Geber_p     = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) t(Geber.f_p(X_,ddt.X,x))))
	E.Geber_p   = t(matrix(apply(Geber_p,2,mean),ncol=nrow(X_)))
	q05.Geber_p = t(matrix(apply(Geber_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.Geber_p = t(matrix(apply(Geber_p,2,quantile,p=0.95),ncol=nrow(X_)))

	## compute relative Geber
	prop.Geber_p     = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) t(prop.Geber.f_p(X_,ddt.X,x))))
	E.prop.Geber_p   = t(matrix(apply(prop.Geber_p,2,mean),ncol=nrow(X_)))
	q05.prop.Geber_p = t(matrix(apply(prop.Geber_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.prop.Geber_p = t(matrix(apply(prop.Geber_p,2,quantile,p=0.95),ncol=nrow(X_)))
	
	## compute mean effects and prop tot Geber
	mean.ddx.Yhat_p    = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) mean.ddx.f_p(X_,x)))
	prop.tot.Geber_p   = t(apply(Omega_p_chain[,-c(1:3)],1,function(x) prop.tot.Geber.f_p(X_,ddt.X,x)))

	## de-scale predictions
	MaP.Yhat_p      = mean_y + sd_y * MaP.Yhat_p
	E.Yhat_p        = mean_y + sd_y * E.Yhat_p
	q05.Yhat_p      = mean_y + sd_y * q05.Yhat_p
	q95.Yhat_p      = mean_y + sd_y * q95.Yhat_p
	#
	MaP.ddx.Yhat_p  = sd_y * t(1/sd_x * t(MaP.ddx.Yhat_p))
	E.ddx.Yhat_p    = sd_y * t(1/sd_x * t(E.ddx.Yhat_p))
	q05.ddx.Yhat_p  = sd_y * t(1/sd_x * t(q05.ddx.Yhat_p))
	q95.ddx.Yhat_p  = sd_y * t(1/sd_x * t(q95.ddx.Yhat_p))
	#
	E.Geber_p       = sd_y * t(1/sd_x * t(E.Geber_p))
	q05.Geber_p     = sd_y * t(1/sd_x * t(q05.Geber_p))
	q95.Geber_p     = sd_y * t(1/sd_x * t(q95.Geber_p))
	#
	mean.ddx.Yhat_p = sd_y/sd_x*mean.ddx.Yhat_p
	
	## store
	results_p         = list("E.Yhat_p"        = E.Yhat_p, 
    						"q05.Yhat_p"       = q05.Yhat_p,    
    						"q95.Yhat_p"       = q95.Yhat_p,   
    						"E.ddx.Yhat_p"     = E.ddx.Yhat_p,
    						"q05.ddx.Yhat_p"   = q05.ddx.Yhat_p,
    						"q95.ddx.Yhat_p"   = q95.ddx.Yhat_p,
    						"E.Geber_p"        = E.Geber_p,   
    						"q05.Geber_p"      = q05.Geber_p,  
    						"q95.Geber_p"      = q95.Geber_p,
    						"E.prop.Geber_p"   = E.prop.Geber_p,   
    						"q05.prop.Geber_p" = q05.prop.Geber_p,  
    						"q95.prop.Geber_p" = q95.prop.Geber_p,
							"mean.ddx.Yhat_p"  = mean.ddx.Yhat_p,
							"prop.tot.Geber_p" = prop.tot.Geber_p 
							)
	return(results_p)
}

## plot.model_p
## goal: 
plot.model_p = function(t,t_,Y,results_p,col="red",xlab=c("Time","Time","Time"),ylab=c("Dynamics","Effects","Contributions"),legend=NULL)
{
	attach(results_p,warn.conflicts=F)
	N = ncol(E.ddx.Yhat_p)
	# par(mfrow=c(3,1),mar=c(2,5,1,1),cex.lab=1.5)
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(Y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=ylab[1])
	points(t_,Y,pch=16)
	polygon(c(t_,rev(t_)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_p,col=col)
	if(!is.null(legend)) legend("topright",legend=legend[1],bty="n",cex=1)
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=ylab[2])
	for(j in 1:N) lines(t_, E.ddx.Yhat_p[,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(legend)) legend("topright",legend=legend[2],bty="n",cex=1)
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=ylab[3])
	for(j in 1:N) lines(t_, E.Geber_p[,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(legend)) legend("topright",legend=legend[3],bty="n",cex=1)
	#
	# par(mfrow=c(1,1))
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
	results_o = list()
	for(i in 1:ncol(Y))
	{
		DFT_        = fit.DFT(Y[,i],t,t_,W_o,log)
		results_o[[i]] = list("MaP.Yhat_o"     = DFT_[,1],
						      "MaP.ddt.Yhat_o" = DFT_[,2],	
						      "E.Yhat_o"       = DFT_[,1],
						      "E.ddt.Yhat_o"   = DFT_[,2],
						      "q05.Yhat_o"     = DFT_[,1],		
	                          "q05.ddt.Yhat_o" = DFT_[,2],
						      "q95.Yhat_o"     = DFT_[,1],		
	                          "q95.ddt.Yhat_o" = DFT_[,2])
	}
	return(results_o)
}

## fit.NODEGM_o
fit.NODEGM_o = function(t,t_,Y,W_o,sd1_o,sd2_o,N_o,log=F,logMar=F)
{
	results_o = list()
	for(i in 1:ncol(Y))
	{
		ensembleList   = fit.model_o(t,Y[,i],W_o,sd1_o,sd2_o,N_o,log,logMar)
		results_o[[i]] = predict.model_o(t_,ensembleList)
	}
	return(results_o)
}

## fit.NODEGM_p
fit.NODEGM_p = function(X,ddt.X,Y,W_p,sd1_p,sd2_p,N_p,log=F,logMar=F)
{
	results_p = list()
	for(i in 1:ncol(Y))
	{
		Omega_p_chain = fit.model_p(X,Y[,i],W_p,sd1_p,sd2_p,N_p,log,logMar)
		results_p[[i]]   = c(summary.model_p(X,ddt.X,Omega_p_chain))
	}
	return(results_p)
}

##
plot.NODEGM_o = function(t,t_,Y,results_o,col="red",xlab=c("","Time"),ylab=c("State","Dynamics"),legend=NULL)
{
	if (is.null(legend)) legend = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
	par(mar=c(4,4,0,0),oma=c(1,1,1,1))
	# par(mar=c(1,1,1,1)*2,oma=c(1,1,1,1))
	layout(matrix(1:(ncol(Y)*2),byrow=F,nrow=2))
	for (i in 1:ncol(Y)) plot.model_o(t,t_,Y[,i],model_o[[i]],col[i],xlab,if (i == 1) ylab else rep("",2),legend[1:2+(i-1)*2])
}

##
plot.NODEGM_p = function(t,t_,Y,results_p,col="red",xlab=c("","","Time"),ylab=c("Dynamics","Effects","Contributions"),legend=NULL)
{
	if (is.null(legend)) legend = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
	par(mar=c(4,4,0,0),oma=c(1,1,1,1))
	# par(mar=c(1,1,1,1)*2,oma=c(1,1,1,1))
	layout(matrix(1:(ncol(Y)*3),byrow=F,nrow=3))
	for (i in 1:ncol(Y)) plot.model_p(t,t_,Y[,i],model_p[[i]],col[i],xlab,if (i == 1) ylab else rep("",3),legend[1:3+(i-1)*3])
}

## diamondPlot.NODEGM_p
## plot mean effects and contributions
diamondPlot.NODEGM_p = function(outputList,labels=NULL)
{
	N = length(outputList)
	par(cex.lab=1,cex.axis=1,cex.main=1,mar=c(1,1,1,1)*2,oma=c(1,5,1,5))
	layout(matrix(1:(N*2),byrow=T,ncol=2))
	mainVect  = matrix(c("a.","b.","c.","d.","e.","f.","g.","h.","i.")[1:(N*2)],byrow=T,ncol=2)
	labels_   = rep("",N)
	#
	for(i in 1:N)
	{
		y = outputList[[i]]$prop.tot.Geber_p
		util.plot.Hdiamond(y,labels=labels,main= if (i == 1) "Total contributions"  else "" )
		legend("topright",legend=mainVect[i,1],bty="n",cex=1)
		#
		y = outputList[[i]]$mean.ddx.Yhat_p
		util.plot.Hdiamond(y,labels=labels_,main= if (i == 1) "Mean effects"  else "" )
		legend("topright",legend=mainVect[i,2],bty="n",cex=1)
	}
	#
	par(mfrow=c(1,1))
}

#
###
