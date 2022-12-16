#########################
## NODE_RMT_functions.r ##
##########################

## goal: 
# - perfrom residual minimisation training on the time series data
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


###############
## UTILITARY ##
###############

## goal: 

## .plot.DIN
#
## goal: plot the dynamical interaction network of the system
#
## input:
## effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
## weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
## labels     - vector - vector of the names of the variables in the matrix
#
## output:
.plot.DIN = function(effectsMat,weightsMat,labels)
{
	## dimensions
	N = dim(effectsMat)[1]

	## scale effects and contributions
	effectsMat = (effectsMat>0)*1

	## angles
	theta = seq(0,2*pi,(2*pi)/N)
	x = cos(theta)
	y = sin(theta)
	x_ = cos(theta+(2*pi*0.05))
	y_ = sin(theta+(2*pi*0.05))
	x__ = 1.25*cos(theta+(2*pi*0.025))
	y__ = 1.25*sin(theta+(2*pi*0.025))

	## plot interactions
	plot(x=c(-1:1)*2,y=c(-1:1)*2,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	for(i in 1:N)
	{
		for(j in 1:N)
		{
			color_ = if(effectsMat[i,j]>0){"green"}else{"red"}
			# points(x__[i],y__[i],cex=30/N)
			text(x__[i],y__[i],labels=labels[i])
			if(weightsMat[i,j]*10>0)
			{
				arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
			}
		}
	}
}

## .plot.diamondplot
#
## goal: plot distributionns of the variables in a table
#
# input:
# df - matrix - matrix containing variables to be plotted by columns
#
## output:
## none
.plot.diamond = function(df,y_lim=NULL,colVect=NULL,y_lab="")
{
	## initiate graphical parameters
	alpha = 0.2 # scaling factor for the plot
	n_col = ncol(df)
	n_row = nrow(df)
	x = 1:n_col
	x_lim = c(min(x)-alpha,max(x)+alpha)
	tmp = as.vector(apply(df,2,function(x)density(x)$x))
	if(is.null(y_lim)) y_lim = c(min(tmp),max(tmp));
	if(is.null(colVect)) colVect = rainbow(n_col,alpha=0.5) 

	## plot
	plot(rep(1,n_row),df[,1],xlim=x_lim,ylim=y_lim,cex=0,xlab="",ylab=y_lab)
	for(j in 1:n_col)
	{
		## coordindates
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

## anchorSampling
#
## goal: function to perform anchor sampling of a given posterior distribution
#
## input:
## n         - int - number of samples to take
## fn_prior  - fn  - function to sample from prior distribution of parameters 
## fn_argmax - fn  - function that optimises posterior distribution given parameters sampled from the priors
#
## output:
## none
anchorSample = function(n,fn_prior,fn_argmax)
{
	ensemble = NULL
	for(i in 1:n)
	{
		message(paste(i,"/",n,sep=""))
		Omega_0 = fn_prior()
		Omega_0 = fn_argmax(Omega_0)
		ensemble = rbind(ensemble,Omega_0)
	}
	return(ensemble)
}

#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## goal: functions for the observation model

## f_o
## goal: interpolate state variable
#
## input:
## t - float - time step in arbitrary units
## Omega - vect - vector of parameters 
#
## output:
## float - value of the interpolation at time t
f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o
## goal: time derivative of the interpolated variable
#
## input:
## t - float - time step in arbitrary units
## Omega - vect - vector of parameters 
#
## output:
## float - value of the derivative wtr to t of the interpolation at time t
ddt.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega_o.f_o
#
## goal: derivative of the interpolation wtr to each network parameter
#
## input:
## t - float - time step in arbitrary units
## Omega - vect - vector of parameters 
#
## output:
## vect - value of the derivative of the interpolation wtr to each parameter in the network
ddOmega_o.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	dfdOmega1 = sin(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
	dfdOmega3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
}

## ddOmega_o.ddt.f_o
#
## goal: derivative of the time derivative of the interpolation wtr to each network parameter
#
## input:
## t - float - time step in arbitrary units
## Omega - vect - vector of parameters 
#
## output:
## vect - value of the derivative of the time derivative of the interpolation wtr to each parameter in the network
ddOmega_o.ddt.f_o = function(t,Omega) 
{	
	Omega = matrit(Omega,ncol=3)
	ddOmega_o1. = pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))
	ddOmega_o2. = Omega[,1]*(pi*cos(pi*(t*Omega[,2] + Omega[,3])) - pi*Omega[,2]*t*sin(pi*(t*Omega[,2] + Omega[,3])))
	ddOmega_o3. = -Omega[,1]*(pi*Omega[,2]*1*sin(pi*(t*Omega[,2] + Omega[,3])))
	return(c(ddOmega_o1.,ddOmega_o2.,ddOmega_o3.))
}

## evaluate functions across multiple time steps
#
## goal: evaluate functions across multiple time steps
#
## input:
## t - float - time step in arbitrary units
## Omega - vect - vector of parameters 
# 
## output:
## values of the corresponding functions evaluated across multiple time steps
f_o.eval           = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval       = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega_o.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega_o.f_o(x,Omega))

# ## logPost_o
# #
# ## goal: evaluate the log posterior of the observation model given the observed data
# #
# ## input:
# ## t       - float - time step in arbitrary units
# ## Y       - vect  - vect containing observations
# ## Omega_o - vect  - vector of parameters 
# #
# ## output:
# ## float   - value of the log of the posterior of observation model
# logPost_o = function(t,Y,Omega_o)
# {
# 	## predict
# 	Ybar_o = f_o.eval(t,Omega_o)
# 
# 	## residuals
# 	res_o = Y - Ybar_o
# 
# 	## Likelihood
# 	logLik_o = - sum(res_o^2)/(sd1_o^2)
# 
# 	## prior
# 	logPrior_o = - sum(Omega_o^2)/(sd2_o^2)
# 
# 	## posterior
# 	logPost_o  = logLik_o + logPrior_o
# 
# 	## terminate
# 	return(logPost_o)
# }
# 
# ddOmega_o.logPost_o = function(t,Y,Omega_o)
# {
# 	## predict
# 	Ybar_o = f_o.eval(t,Omega_o)
# 
# 	## derivative of residuals
# 	res_o = (Y - Ybar_o)
# 	ddOmega_o.res_o = - ddOmega_o.f_o.eval(t,Omega_o)
# 
# 	## derivative of likelihood
# 	ddOmega_o.logLik_o = - 2 * ddOmega_o.res_o%*%res_o/(sd1_o^2)
# 
# 	## derivative of prior
# 	ddOmega_o.logPrior_o = - 2 * Omega_o/(sd2_o^2)
# 
# 	## derivative of posterior
# 	ddOmega_o.logPost_o = ddOmega_o.logLik_o + ddOmega_o.logPrior_o
# 
# 	## terminate
# 	return(ddOmega_o.logPost_o)
# }

## logPost_o
## goal: compute the log marginal posterior density of a parameter vector
#
## input:
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters 
#
## output:
## float   - value of the log marginal posterior density
logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## residuals
	res_o = Y - Ybar_o

	## Likelihood
	logLik_o = - 0.5 * length(Y) * log(0.5 * sum(res_o^2) + 1)

	## prior
	logPrior_o = - 0.5 * length(Omega_o) * log(0.5 * sum(Omega_o^2) + 1)

	## posterior
	logPost_o  = logLik_o + logPrior_o

	## terminate
	return(logPost_o)
}

## ddOmega_o.logPost_o
#
## goal: compute the derivate of the log marginal posterior density wtr to each parameter
#
## input:
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - vector of parameters 
#
## output:
## vector  - value of the derivative for eah parameter
ddOmega_o.logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## derivative of residuals
	res_o = (Y - Ybar_o)
	ddOmega_o.res_o = -ddOmega_o.f_o.eval(t,Omega_o)

	## derivative of likelihood
	ddOmega_o.logLik_o = - 0.5 * length(Y) * 1/(0.5 * sum(res_o^2) + 1) * 0.5 * ddOmega_o.res_o%*%res_o

	## derivative of prior
	ddOmega_o.logPrior_o = - 0.5 * length(Omega_o) * 1/(0.5 * sum(Omega_o^2) + 1) * Omega_o

	## derivative of posterior
	ddOmega_o.logPost_o = ddOmega_o.logLik_o + ddOmega_o.logPrior_o

	## terminate
	return(ddOmega_o.logPost_o)
}

## argmax.logPost_o
#
# goal: maximise the log posterior density of the parameters
#
## input:
## t       - float - time step in arbitrary units
## Y       - vect  - vect containing observations
## Omega_o - vect  - initial vector of parameters 
#
## output:
## Omega_o - vect  - vector of parameters that maximise (locally) the posterior
argmax.logPost_o = function(t,Y,Omega_o)
{
	error_     = function(x) -logPost_o(t,Y,x)
	graderror_ = function(x) -ddOmega_o.logPost_o(t,Y,x)
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

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## activation functions
f_sigma     = function(x) 1/(1+exp(-x)) 
ddu.f_sigma = function(x) f_sigma(x) * (1 - f_sigma(x))
# f_sigma     = function(x) x
# ddu.f_sigma = function(x) 1
# f_sigma     = function(x) (x>0)*x
# ddu.f_sigma = function(x) (x>0)
# f_sigma     = function(x) exp(-x^2)
# ddu.f_sigma = function(x) -2*x*exp(-x^2)

## ddt.f_p 
#
# goal: process model defined as a neural ordinary differential equation (NODE) 
#       to approximate dynamics of each state variable
#
## input:
## x     - vector - vector of input variables
## Omega - vector - vector of parameters of the NODE
#
## output:
## float - value of the NODE
ddt.f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,ncol=2 + length(x))
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
}

## ddOmega_p.ddt.f_p
#
## goal: derivative of the process model wtr to parameters
#
## input:
## x      - vector - vector of input variables
## Omega  - vector - vector of parameters of the NODE
#
## output:
## vector - value of the derivative of the NODE wtr to each parameter
ddOmega_p.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	ddOmega_p1. = f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p2. = Omega[,1]*1*ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p3. = Omega[,1]*rep(1,W_p)%*%t(x)*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
	return(c(ddOmega_p1.,ddOmega_p2.,ddOmega_p3.))
}

## ddx.ddt.f_p
#
## goal: derivative of the process model wtr to state variables
#
## input:
## x      - vector - vector of input variables
## Omega  - vector - vector of parameters of the NODE
#
## output:
## vector - value of the derivative of the NODE wtr to each state variable
ddx.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,ncol=2 + length(x))
	ddx. = Omega[,1]%*%(Omega[,-c(1:2)]*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)))
	return(ddx.)
}

## evaluate process functions across several timesteps
#
## input:
## X      - matrix - matrix of variables across all time steps
## Omega  - vector - vector of parameters of the network 
#
## output:
## vector - value of the function evaluated across all time steps
ddt.f_p.eval           = function(X,Omega) apply(X,1,function(x) ddt.f_p(x,Omega))
ddx.ddt.f_p.eval       = function(X,Omega) apply(X,1,function(x) ddx.ddt.f_p(x,Omega))
ddOmega_p.ddt.f_p.eval = function(X,Omega) apply(X,1,function(x) ddOmega_p.ddt.f_p(x,Omega))

## logPost_p
#
# goal: log posterior density of the process model (first level of inference)
#
## input:
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
#
## output:
## float      - value of the posterior density givent the interpolated data and parameter vector
logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
  ## predict
  ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)

  ## residuals
  res_p = ddt.Ybar_o - ddt.Ybar_p

  ## Likelihood
  logLik_p   = - sum(res_p^2)/(sd1_p^2)

  ## prior
  logPrior_p = - sum(Omega_p^2)/(sd2_p^2)

  ## posterior
  logPost_p  = logLik_p + logPrior_p

  ## terminate
  return(logPost_p)
}

## ddOmega_p.logPost_p
#
# goal: derivative of the log posterior density of the process model wtr to parameters
#
## input:
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
#
## output:
## vector     - value of the derivative of the posterior density wtr to each parameter 
ddOmega_p.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
  ## predict
  ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)

  ## derivative of residuals
  res_p           = ddt.Ybar_o - ddt.Ybar_p
  ddOmega_p.res_p = - ddOmega_p.ddt.f_p.eval(Ybar_o,Omega_p)

  ## derivative of likelihood
  ddOmega_p.logLik_p = - 2 * ddOmega_p.res_p%*%res_p/(sd1_p^2)

  ## derivative of prior
  ddOmega_p.logPrior_p = - 2 * Omega_p /(sd1_p^2)

  ## derivative of posterior
  ddOmega_p.logPost_p = ddOmega_p.logLik_p + ddOmega_p.logPrior_p

  ## terminate
  return(ddOmega_p.logPost_p)
}

## argmax.logPost_p
#
# goal: find the maximum of the log posterior density of the process model
#
## input:
## Ybar_o     - matrix - matrix containing interpolated time series of the variables 
## ddt.Ybar_o - matrix - matrix containing time derivative of interpolated variables
## Omega_p    - vector - parameter vector 
#
## output:
## vector     - vector of parameters that maximise locally the posterior density
argmax.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
	error_     = function(x) -logPost_p(Ybar_o,ddt.Ybar_o,x)
	graderror_ = function(x) -ddOmega_p.logPost_p(Ybar_o,ddt.Ybar_o,x)
	Omega_p  = optim(par=Omega_p,
				fn=error_,
				gr=graderror_,
				method="BFGS"#,
				# control=list("trace"=1,"REPORT"=1,"maxit"=100)
				)$par
	return(Omega_p)
}

# ## compute equilibrium densities
# argmin.ddt.f_p = function(Y_0,Omega)
# {
#   ## minimise ddt.f_p
#   Omega = matrix(Omega,ncol=N,byrow=F)
#   fn = function(x) ddt.f_p(x, Omega[,1])^2 + ddt.f_p(x, Omega[,2])^2 + ddt.f_p(x, Omega[,3])^2
#   res = optim(par=Y_0, fn=fn)
#   return(c(res$value,res$par))
# }
# #
# tmp = apply(ensemble,1, function(x) argmin.ddt.f_p(rep(0,N),x))
# tmp = tmp[,which(tmp[1,]<=0.001)][-1,]
# for(i in 1:N) 
# {
#    	torm = which(abs(tmp[i,])>=3) 
# 	if( length(torm) != 0 )
# 	{
# 		tmp = tmp[,-torm]
# 	}
# }
# argmin.ddt.f_p_q0.05 = matrix(apply(tmp,1,quantile,probs=0.05),nrow=N,byrow=T)
# argmin.ddt.f_p_q0.50 = matrix(apply(tmp,1,quantile,probs=0.50),nrow=N,byrow=T)
# argmin.ddt.f_p_q0.95 = matrix(apply(tmp,1,quantile,probs=0.95),nrow=N,byrow=T)
# #
# tmp = t(tmp); colnames(tmp) = colnames(Y)
# boxplot(tmp,ylab="Density")

#
###
