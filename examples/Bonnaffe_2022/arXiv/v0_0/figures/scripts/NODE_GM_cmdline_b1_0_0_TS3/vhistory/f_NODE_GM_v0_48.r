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

## next:
##            - fix number of param issue in fit.model_p
##            - add log of response variable in process model
##            - make sure description of functions are accurate
##            - review name of functions
##            - allow for different values of priors for each time series
##            - fix number of parameters in model
##            - allow to switch activation functions

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
# log    - bool   - whether to log the response variable (e.g. if its support is positive)
# logMar - bool   - whether to optimise the log posterior or marginal
# 
## next:
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
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
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

## fit.NODEGM_o ##
## goal: compute interpolation/extrapolation and derivative of response variable of observation model using Single Layer Perceptron (SLPs)
##       procedure is applied to all response variables in matrix
## out:  returns list containing for each response variable the scaling parameters and 
##       the inferred interpolation/extrapolation and temporal derivatives
# t       - vector - time steps
# t_      - vector - interpolated/extrapolated time steps
# Y       - matrix - values of response variables
# W       - int    - number of periodic elements in the series
# sd_1    - float  - standard deviation of likelihood
# sd_2    - vector - standard deviation of prior distributions
# N_e     - int    - number of elements to sample
# log     - bool   - whether to log the response variable (e.g. if positive support only)
# logMar  - bool   - whether to use the log marginal instead of the log posterior
fit.NODEGM_o = function(t,t_,Y,W,sd_1,sd_2,N_e,log=F,logMar=F)
{
	results = list()
	for(i in 1:ncol(Y))
	{
		ensembleList = fit.model_o(t,Y[,i],W,sd_1,sd_2,N_e,log,logMar)
		results[[i]] = predict.model_o(t_,ensembleList)
	}
	return(results)
}

## plot.NODEGM_o ##
## goal: visualise results of fitting  
# t         - vector - time steps
# t_        - vector - interpolated/extrapolated time steps
# Y         - vector - values of response variables
# results_o - list   - predictions of observation model (i.e. interpolated/extrapolated response and temporal derivatives)
# ...
plot.NODEGM_o = function(t,t_,Y,results,col="red",xlab=c("","Time"),ylab=c("State","Dynamics"),index=NULL)
{
	## initiate
	if (is.null(index)) index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
	par(mar=c(4,4,0,0),oma=c(1,1,1,1),cex.lab=1.25)
	layout(matrix(1:(ncol(Y)*2),byrow=F,nrow=2))

	## plot
	for (i in 1:ncol(Y)) 
	{
		plot.model_o(t,t_,Y[,i],results[[i]],col[i],xlab,if (i == 1) ylab else rep("",2),index[1:2+(i-1)*2])
	}
	#
	par(mfrow=c(1,1),mar=c(4,4,3,3))
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
# f_sigma     = function(x) exp(x)
# ddu.f_sigma = function(x) exp(x)
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
f_sigma_p     = sigmo
ddx.f_sigma_p = ddx.sigmo

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
logLik_p = function(X,Y,Omega,sd_1,sd_2)
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
# log    - bool   - whether to log the response variable (e.g. if its support is positive)
# logMar - bool   - whether to optimise the log posterior or marginal
fit.model_p = function(X,Y,W,sd_1,sd_2,N_e,log=F,logMar=F)
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
	ensemble= NULL 
	for(k in 1:N_e)    
	{   
	    Omega_0 = rnorm(W*(2+ncol(X_))*2,0,sd_2)
		if(logMar==F)
		{
	    	Omega_f   = argmax.logPost_p(X_,Y_,Omega_0,sd_1,sd_2)
	    	logPost_0 =        logPost_p(X_,Y_,Omega_0,sd_1,sd_2)
	    	logPost_f =        logPost_p(X_,Y_,Omega_f,sd_1,sd_2)
		} else
		{
			Omega_f   = argmax.logMar_p(X_,Y_,Omega_0)
	    	logPost_0 =        logMar_p(X_,Y_,Omega_0)
	    	logPost_f =        logMar_p(X_,Y_,Omega_f)
		}
	    ensemble = rbind(ensemble,c(logPost_f,Omega_f))
	    message(paste(k,"/",N_e,"\t",
	            format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
	            format(round(logPost_f,2),nsmall=2),sep=""))
		message(paste("loglik:\t\t",format(round(logLik_p(X_,Y_,Omega_f,sd_1)),nsmall=2)))
		message(paste("logprior:\t",format(round(logPrior_p(Omega_f,sd_2))    ,nsmall=2)))
		message(paste("r2:\t\t",    format(round(r2_p(X_,Y_,Omega_f),2)         ,nsmall=2)))

	}
	
	## terminate
	message("") 
	message(paste("MaP: ",max(ensemble[,1]),sep=""))
	message("") 
	results = list("metaData" = list("mean_x" = mean_x,
								     "sd_x"   = sd_x,
									 "mean_y" = mean_y,
								     "sd_y"   = sd_y,
									 "log"    = log),
				   "ensemble" = ensemble
				  )
	return(results)
}

## predict.model_p ##
## goal: predict response for values of X explanatory variables
# X             - matrix - explanatory variables (e.g. state variables)
# Y             - vector - response variable
# ensembleLisst - list   - list containing metadata to scale/de-scale variables and ensemble of parameter vectors
predict.model_p = function(X,ensembleList)
{
	## unpack results
	metaData = ensembleList$metaData
	ensemble = ensembleList$ensemble

	## standardisation parameters
	mean_x = metaData[["mean_x"]]
	sd_x   = metaData[["sd_x"]]
	mean_y = metaData[["mean_y"]]
	sd_y   = metaData[["sd_y"]]
	log    = metaData[["log"]]

	## standardise predictive variables
	mean_x = apply(X,2,mean)
	sd_x   = apply(X,2,sd)
	X_     = t((t(X)-mean_x)/sd_x)
	
	## maximum a posteriori
	idx            = which.max(ensemble[,1])
	MaP.Omega_p    = ensemble[idx,-1]
	MaP.Yhat_p     = f_p.eval(X_,MaP.Omega_p)
	MaP.ddx.Yhat_p = t(ddx.f_p.eval(X_,MaP.Omega_p)) 
	
	## compute predicted response 
	Yhat_p         = t(apply(ensemble[,-1],1,function(x) f_p.eval(X_,x)))
	E.Yhat_p       = apply(Yhat_p,2,mean)
	q05.Yhat_p     = apply(Yhat_p,2,quantile,p=0.05)
	q95.Yhat_p     = apply(Yhat_p,2,quantile,p=0.95)

	## compute derivative of predicted response wtr explanatory variables
	ddx.Yhat_p     = t(apply(ensemble[,-1],1,function(x) ddx.f_p.eval(X_,x)))
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
	results_p  = list("E.Yhat_p"       = E.Yhat_p, 
    		 	      "q05.Yhat_p"     = q05.Yhat_p,    
    		 	      "q95.Yhat_p"     = q95.Yhat_p,   
    		 	      "E.ddx.Yhat_p"   = E.ddx.Yhat_p,
    		 	      "q05.ddx.Yhat_p" = q05.ddx.Yhat_p,
    		 	      "q95.ddx.Yhat_p" = q95.ddx.Yhat_p)
    return(results_p)
}

## summary.model_p ##
## goal: predict response variable for values of X explanatory variables, their dynamics, contributions (i.e. Geber)
# X             - matrix - explanatory variables (e.g. state variables)
# Y             - vector - response variable
# ensembleLisst - list   - list containing metadata to scale/de-scale variables and ensemble of parameter vectors
summary.model_p = function(X,ddt.X=NULL,ensembleList)
{
	## unpack results
	metaData = ensembleList$metaData
	ensemble = ensembleList$ensemble

	## standardisation parameters
	mean_x = metaData[["mean_x"]]
	sd_x   = metaData[["sd_x"]]
	mean_y = metaData[["mean_y"]]
	sd_y   = metaData[["sd_y"]]
	log    = metaData[["log"]]

	## standardise predictive variables
	X_     = t((t(X)-mean_x)/sd_x)

	## maximum a posteriori
	idx              = which.max(ensemble[,1])
	MaP.Omega_p      = ensemble[idx,-1]
	MaP.Yhat_p       = f_p.eval(X_,MaP.Omega_p)
	MaP.ddx.Yhat_p   = t(ddx.f_p.eval(X_,MaP.Omega_p))
	
	## compute predicted response 
	Yhat_p           = t(apply(ensemble[,-1],1,function(x) f_p.eval(X_,x)))
	E.Yhat_p         = apply(Yhat_p,2,mean)
	q05.Yhat_p       = apply(Yhat_p,2,quantile,p=0.05)
	q95.Yhat_p       = apply(Yhat_p,2,quantile,p=0.95)
	
	## compute derivative of predicted response wtr explanatory variables
	ddx.Yhat_p       = t(apply(ensemble[,-1],1,function(x) ddx.f_p.eval(X_,x)))
	E.ddx.Yhat_p     = t(matrix(apply(ddx.Yhat_p,2,mean),ncol=nrow(X_)))
	q05.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p,2,quantile,p=0.95),ncol=nrow(X_)))

	## compute Geber
	if (is.null(ddt.X)) ddt.X = Yhat_p 
	Geber_p          = t(apply(ensemble[,-1],1,function(x) t(Geber.f_p(X_,ddt.X,x))))
	E.Geber_p        = t(matrix(apply(Geber_p,2,mean),ncol=nrow(X_)))
	q05.Geber_p      = t(matrix(apply(Geber_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.Geber_p      = t(matrix(apply(Geber_p,2,quantile,p=0.95),ncol=nrow(X_)))

	## compute relative Geber
	prop.Geber_p     = t(apply(ensemble[,-1],1,function(x) t(prop.Geber.f_p(X_,ddt.X,x))))
	E.prop.Geber_p   = t(matrix(apply(prop.Geber_p,2,mean),ncol=nrow(X_)))
	q05.prop.Geber_p = t(matrix(apply(prop.Geber_p,2,quantile,p=0.05),ncol=nrow(X_)))
	q95.prop.Geber_p = t(matrix(apply(prop.Geber_p,2,quantile,p=0.95),ncol=nrow(X_)))
	
	## compute mean effects and prop tot Geber
	mean.ddx.Yhat_p  = t(apply(ensemble[,-1],1,function(x) mean.ddx.f_p(X_,x)))
	prop.tot.Geber_p = t(apply(ensemble[,-1],1,function(x) prop.tot.Geber.f_p(X_,ddt.X,x)))

	## de-scale predictions
	MaP.Yhat_p       = mean_y + sd_y * MaP.Yhat_p
	E.Yhat_p         = mean_y + sd_y * E.Yhat_p
	q05.Yhat_p       = mean_y + sd_y * q05.Yhat_p
	q95.Yhat_p       = mean_y + sd_y * q95.Yhat_p
	#
	MaP.ddx.Yhat_p   = sd_y * t(1/sd_x * t(MaP.ddx.Yhat_p))
	E.ddx.Yhat_p     = sd_y * t(1/sd_x * t(E.ddx.Yhat_p))
	q05.ddx.Yhat_p   = sd_y * t(1/sd_x * t(q05.ddx.Yhat_p))
	q95.ddx.Yhat_p   = sd_y * t(1/sd_x * t(q95.ddx.Yhat_p))
	#
	E.Geber_p        = sd_y * t(1/sd_x * t(E.Geber_p))
	q05.Geber_p      = sd_y * t(1/sd_x * t(q05.Geber_p))
	q95.Geber_p      = sd_y * t(1/sd_x * t(q95.Geber_p))
	#
	mean.ddx.Yhat_p  = sd_y/sd_x*mean.ddx.Yhat_p
	
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
							"prop.tot.Geber_p" = prop.tot.Geber_p,
							"ensemble"         = ensembleList
							)
	return(results_p)
}

## plot.model_p ##
## goal: visualise results of process model 
# t_     - vector - time steps
# Y      - vector - response variable that the process model has been fitted to
# ...
# index  - vector - index of graphs
# legend - vector - legend of different effects and contributions 
plot.model_p = function(t_,Y,results_p,col="red",xlab=c("Time","Time","Time"),ylab=c("Dynamics","Effects","Contributions"),index=NULL,legend=NULL)
{
	## initiate
	attach(results_p,warn.conflicts=F)
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

## fit.NODEGM_p ##
## goal: compute predictions of process model for all response variables in table
# X      - matrix - value of explanatory variables at each time step
# ddt.X  - matrix - temporal dynamics of response variables at each time step
# Y      - matrix - value of response variables at each time step
# sd_1   - int    - standard deviation of likelihood
# sd_2   - vector - standard deviation of prior
# N_e    - int    - number of elements to sample
# log    - bool   - whether to log explanatory variables (to facilitated fitting)
# logMar - bool   - whether to use the log marginal instead of the log posterior
fit.NODEGM_p = function(X,ddt.X,Y,W_p,sd1_p,sd2_p,N_p,log=F,logMar=F)
{
	results = list()
	for(i in 1:ncol(Y))
	{
		ensembleList   = fit.model_p(X,Y[,i],W_p,sd1_p,sd2_p,N_p,log,logMar)
		results[[i]] = c(summary.model_p(X,ddt.X,ensembleList))
	}
	return(results)
}

## plot.NODEGM_p ##
## goal: plot predictions of process model
# t_      - vector - interpolated/extrapolated time steps
# Y       - matrix - value of response variables at each time step
# results - list   - predictions of process model (i.e. predicted response, sensitivity of response to each explanatory variable, ...)
# ...
plot.NODEGM_p = function(t_,Y,results,col="red",xlab=c("","","Time"),ylab=c("Dynamics","Effects","Contributions"),index=NULL,legend=NULL)
{
	## initiate
	par(mar=c(4,4,0,0),oma=c(1,1,1,1),cex.lab=1.25)
	layout(matrix(1:(ncol(Y)*3),byrow=F,nrow=3))
	if (is.null(index)) index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")

	## plot
	for (i in 1:ncol(Y)) 
	{
		if (is.null(legend)) legend_ = paste(colnames(Y),"->",colnames(Y)[i])
		plot.model_p(t_,Y[,i],results[[i]],col[i],xlab,if (i == 1) ylab else rep("",3),index[1:3+(i-1)*3],legend=legend_)
	}
	#
	par(mfrow=c(1,1),mar=c(4,4,3,3))
}

## diamondPlot.NODEGM_p ##
## goal: plot mean effects and contributions
# results - list   - predictions of process model (i.e. predicted response, sensitivity of response to each explanatory variable, ...)
diamondPlot.NODEGM_p = function(results,labels=NULL)
{
	## initiate
	N = length(results)
	par(cex.lab=1,cex.axis=1,cex.main=1,mar=c(1,1,1,1)*2,oma=c(1,5,1,5),cex.lab=1.25)
	layout(matrix(1:(N*2),byrow=T,ncol=2))
	mainVect  = matrix(c("a.","b.","c.","d.","e.","f.","g.","h.","i.")[1:(N*2)],byrow=T,ncol=2)
	labels_   = rep("",N)
	#
	## plot
	for(i in 1:N)
	{
		y = results[[i]]$prop.tot.Geber_p
		util.plot.Hdiamond(y,labels=labels,main= if (i == 1) "Total contributions"  else "" )
		legend("topright",legend=mainVect[i,1],bty="n",cex=1)
		#
		y = results[[i]]$mean.ddx.Yhat_p
		util.plot.Hdiamond(y,labels=labels_,main= if (i == 1) "Mean effects"  else "" )
		legend("topright",legend=mainVect[i,2],bty="n",cex=1)
	}
	#
	par(mfrow=c(1,1),mar=c(4,4,3,3))
}

#
###
