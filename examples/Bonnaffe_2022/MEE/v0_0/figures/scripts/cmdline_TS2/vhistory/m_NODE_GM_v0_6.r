############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

## update
# 10-03-2021 - created version 0.0
# 11-03-2021 - created version 0.1
#            - implemented figures for the paper 
# 29-03-2021 - created version 0.2
#            - combined figures for the paper
# 25-05-2021 - created version 0.3
#            - modularised code to work with any single time series file
# 18-11-2021 - crated version 0_5
#            - cleaned code
# 22-11-2021 - cleaned code
# 22-11-2021 - created v0_6

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## load NODE functions
source("f_NODE_GM_v0_28.r")

## load data
TS = read.table("TS.csv",sep=";",header=T)

## transform positive variables
TS[,c("Hare","Lynx")] = log(TS[,c("Hare","Lynx")])

## plot time series
plot(TS[,1],TS[,2],xlab="Time (days)",ylab="Density (normalised)")
for(j in 1:2)
{
	points(TS[,1],TS[,j+1])
	lines(TS[,1],TS[,j+1],col=rainbow(ncol(TS)-1)[j])
}

## parameters
N_e    = 30     # number of ensemble elements
W_o    = 100     # number of nodes in observation model
W_p    = 20     # number of nodes in process model
sd1_o  = 0.01   # sd in likelihood of residuals of observation model
sd2_o  = 0.1    # sd in priors of observation model
sd1_p  = 0.01   # sd in likelihood of residuals of process model
sd2_p  = 0.1    # sd in priors of process model
alpha  = 50     # scaling factor for observation model
method = 1      # 1 for normal posterior / 2 for marginal posterior

## initiate NODE
n     = nrow(TS)                                  # number of time steps
n_    = n*1                                       # number of interpolated time steps
N     = ncol(TS) - 1                              # number of variables
t     = TS[,1]                                    # vector of time steps
Y     = TS[,-1]                                   # matrix containing time series
t     = (t-min(t))/(max(t)-min(t))*alpha          # normalise time steps
t_    = seq(min(t),max(t),(max(t)-min(t))/n_)     # normalise interpolated time steps
n_    = length(t_)                                # double check number of interpolated time steps
mean_ = apply(Y,2,mean)                           # save mean of each variable in time series
sd_   = apply(Y,2,sd)                             # save sd of each variable in time series
Y     = t(apply(Y,1,function(x)(x-mean_)/sd_))    # standardise each variable

#
###

#######################
## OBSERVATION MODEL ##
#######################

## goal: interpolate the time series

## fit observation model
ensembleList_o = list()
for(i in 1:N)
{
	message(paste("Interpolating",colnames(Y)[i]))
  	## get ensemble (anchor sampling)
	fn_prior  = function() rnorm(3*W_o,0,sd2_o)
	fn_argmax = function(x) argmax.logPost_o(t,Y[,i],x)
	ensembleList_o[[i]] = anchorSample(N_e,fn_prior,fn_argmax)
}

## get expectations 
outputList_o = list()
for(i in 1:N)
{
	## load interpolation ensemble
	ensemble = ensembleList_o[[i]]

	## get Ybar
	f_o.eval_  = apply(ensemble,1,function(x)f_o.eval(t_,x))
	q05.Ybar_o = apply(f_o.eval_,1,quantile,probs=0.05)
	q95.Ybar_o = apply(f_o.eval_,1,quantile,probs=0.95)
	E.Ybar_o   = apply(f_o.eval_,1,mean)
	# Ybar_o_ensemble = t(f_o.eval_)
	
	## get ddt.Ybar
	ddt.f_o.eval_ = apply(ensemble,1,function(x)ddt.f_o.eval(t_,x))
	q05.ddt.Ybar_o = apply(ddt.f_o.eval_,1,quantile,probs=0.05)
	q95.ddt.Ybar_o = apply(ddt.f_o.eval_,1,quantile,probs=0.95)
	E.ddt.Ybar_o  = apply(ddt.f_o.eval_,1,mean)
	# ddt.Ybar_o_ensemble = t(ddt.f_o.eval_)

	## store
	outputList_o[[i]] = list(
		"E.Ybar_o"       = E.Ybar_o,
		"q05.Ybar_o"     = q05.Ybar_o,
		"q95.Ybar_o"     = q95.Ybar_o,
		"E.ddt.Ybar_o"   = E.ddt.Ybar_o,
		"q05.ddt.Ybar_o" = q05.ddt.Ybar_o,
		"q95.ddt.Ybar_o" = q95.ddt.Ybar_o
	)
}

## visualise fit
for(i in 1:N)
{
	## open outputList
	# attach(outputList_o[[i]])
	E.Ybar_o       = outputList_o[[i]][["E.Ybar_o"]]
	q05.Ybar_o     = outputList_o[[i]][["q05.Ybar_o"]]
	q95.Ybar_o     = outputList_o[[i]][["q95.Ybar_o"]]
	E.ddt.Ybar_o   = outputList_o[[i]][["E.ddt.Ybar_o"]]
	q05.ddt.Ybar_o = outputList_o[[i]][["q05.ddt.Ybar_o"]]
	q95.ddt.Ybar_o = outputList_o[[i]][["q95.ddt.Ybar_o"]]
	
	## plot
	par(mfrow=c(2,1),mar=c(4,5,1,4),cex.lab=1)
	#
	## visualise dynamics
	plot(t,Y[,i],type="b",ylim=c(-3,3),xlab="Time",ylab="State")
	polygon(c(t_,rev(t_)),c(q05.Ybar_o,rev(q95.Ybar_o)),col=rainbow(N,alpha=0.2)[i],border=NA)
	lines(t_,E.Ybar_o,col=rainbow(N)[i],lty=1)
	#	
	## visualise derivative
	plot(t_,E.ddt.Ybar_o,ylim=c(-3,3),cex=0,xlab="Time",ylab="Dynamics")
	polygon(c(t_,rev(t_)),c(q05.ddt.Ybar_o,rev(q95.ddt.Ybar_o)),col=rainbow(N,alpha=0.2)[i],border=NA)
	lines(t_,E.ddt.Ybar_o,col=rainbow(N)[i],lty=1)
	#	
	par(mfrow=c(1,1))
}

#
###

###################
## PROCESS MODEL ##
###################

## goal: fit the process model

## format input
Ybar_o = NULL
ddt.Ybar_o = NULL
for (i in 1:N)
{
	Ybar_o = cbind(Ybar_o,outputList_o[[i]][["E.Ybar_o"]])
	ddt.Ybar_o = cbind(ddt.Ybar_o,outputList_o[[i]][["E.ddt.Ybar_o"]])
}

## get ensemble
ensembleList_p = list()
for (i in 1:N)
{
	message(paste("Fitting",colnames(Y)[i]))
	## get ensemble (anchor sampling)
	fn_prior  = function()  rnorm((2+N)*W_p,0,sd2_p)
	fn_argmax = function(x) argmax.logPost_p(Ybar_o,ddt.Ybar_o[,i],x)
	ensemble = anchorSample(N_e,fn_prior,fn_argmax)
	ensembleList_p[[i]] = ensemble
}

## get expectations 
outputList_p = list()
for(i in 1:N)
{
	## get ensemble for a specific time series
	ensemble = ensembleList_p[[i]]
	
	## compute ddt.Ybar_p
	tmp = apply(ensemble,1, function(x) ddt.f_p.eval(Ybar_o,x))
	q05.ddt.Ybar_p = apply(tmp,1,quantile,probs=0.05)
	q95.ddt.Ybar_p = apply(tmp,1,quantile,probs=0.95)
	E.ddt.Ybar_p = apply(tmp,1,mean)
	# ddt.Ybar_p_ensemble = t(tmp)
	
	## compute effects
	tmp = apply(ensemble,1, function(x) ddx.ddt.f_p.eval(Ybar_o,x))
	q05.ddx.ddt.Ybar_p = matrix(apply(tmp,1,quantile,probs=0.05),nrow=n_,byrow=T)
	q95.ddx.ddt.Ybar_p = matrix(apply(tmp,1,quantile,probs=0.95),nrow=n_,byrow=T)
	E.ddx.ddt.Ybar_p = matrix(apply(tmp,1,mean),nrow=n_,byrow=T)
	# ddx.ddt.Ybar_p_ensemble = t(tmp)
	
	## compute mean effects
	# tmp = apply(ensemble, 1, function(x) apply(ddx.ddt.f_p.eval(Ybar_o,x), 1, mean))
	E.mean.ddx.ddt.Ybar_p = apply(E.ddx.ddt.Ybar_p,2,mean)
	# mean.ddx.ddt.Ybar_p_ensemble = t(tmp)
	
	## compute contributions
	tmp = apply(ensemble, 1, function(x) ddx.ddt.f_p.eval(Ybar_o,x) * t(ddt.Ybar_o))
	q05.C = matrix(apply(tmp,1,quantile,probs=0.05),nrow=n_,byrow=T)
	q95.C = matrix(apply(tmp,1,quantile,probs=0.95),nrow=n_,byrow=T)
	E.C = matrix(apply(tmp,1,mean),nrow=n_,byrow=T)
	# C_ensemble = t(tmp)
	
	## compute relative contributions
	tmp = apply(ensemble, 1, function(x) apply(ddx.ddt.f_p.eval(Ybar_o,x) * t(ddt.Ybar_o), 2, function(x) x^2/sum(x^2)))
	q05.propC = matrix(apply(tmp,1,quantile,probs=0.05,na.rm=T),nrow=n_,byrow=T)
	q95.propC = matrix(apply(tmp,1,quantile,probs=0.95,na.rm=T),nrow=n_,byrow=T)
	E.propC = matrix(apply(tmp,1,mean),nrow=n_,byrow=T)
	# propC_ensemble = t(tmp)
	
	## compute relative total contributions
	tmp = apply(ensemble, 1, function(x) apply(ddx.ddt.f_p.eval(Ybar_o,x) * t(ddt.Ybar_o), 1, function(x) sum(x^2)))
	tmp = apply(tmp, 2, function(x) x/sum(x)) # normalise
	q05.prop.totC = apply(tmp,1,quantile,probs=0.05)
	q95.prop.totC = apply(tmp,1,quantile,probs=0.95)
	E.prop.totC = apply(tmp,1,mean)
	# prop.totC_ensemble = t(tmp)
	  
	## outputs
	outputList_p[[i]] = list(
		"E.ddt.Ybar_p"       = E.ddt.Ybar_p,
		"q05.ddt.Ybar_p"     = q05.ddt.Ybar_p,
		"q95.ddt.Ybar_p"     = q95.ddt.Ybar_p,
		"E.ddx.ddt.Ybar_p"   = E.ddx.ddt.Ybar_p,
		"q05.ddx.ddt.Ybar_p" = q05.ddx.ddt.Ybar_p,
		"q95.ddx.ddt.Ybar_p" = q95.ddx.ddt.Ybar_p,
		"E.C"                = E.C,
		"q05.C"              = q05.C,
		"q95.C"              = q95.C,
		"E.propC"            = E.propC,
		"q05.propC"          = q05.propC,
		"q95.propC"          = q95.propC
	) 
}

## visualise results
for(i in 1:N)
{
	## open outpuList
	# attach(outputList_o[[i]])
	# attach(outputList_p[[i]])
	#
	E.Ybar_o           = outputList_o[[i]][["E.Ybar_o"]]
	E.ddt.Ybar_o       = outputList_o[[i]][["E.ddt.Ybar_o"]]
	E.ddt.Ybar_p       = outputList_p[[i]][["E.ddt.Ybar_p"]]
	q05.ddt.Ybar_p     = outputList_p[[i]][["q05.ddt.Ybar_p"]]
	q95.ddt.Ybar_p     = outputList_p[[i]][["q95.ddt.Ybar_p"]]
	E.ddx.ddt.Ybar_p   = outputList_p[[i]][["E.ddx.ddt.Ybar_p"]]
	q05.ddx.ddt.Ybar_p = outputList_p[[i]][["q05.ddx.ddt.Ybar_p"]]
	q95.ddx.ddt.Ybar_p = outputList_p[[i]][["q95.ddx.ddt.Ybar_p"]]
	E.C                = outputList_p[[i]][["E.C"]]
	q05.C              = outputList_p[[i]][["q05.C"]] 
	q95.C              = outputList_p[[i]][["q95.C"]]
	E.propC            = outputList_p[[i]][["E.propC"]]
	q05.propC          = outputList_p[[i]][["q05.propC"]]
	q95.propC          = outputList_p[[i]][["q95.propC"]]
	
	## plot dynamics
	par(mfrow=c(3,1),mar=c(4,5,1,4),cex.lab=1)
	#
	## plot dynamics
	plot(t_,E.ddt.Ybar_o,ylim=c(-3,3),type="l",xlab="Time",ylab="Dynamics")
	polygon(c(t_,rev(t_)),c(q05.ddt.Ybar_p,rev(q95.ddt.Ybar_p)),border=NA,col=rainbow(N,alpha=0.2)[i])
	lines(t_,E.ddt.Ybar_p,type="l",col=rainbow(N)[i])
	lines(t_,rep(0,n_),lty=2)
	#
	#
	## plot effects 
	X_q0.05 = q05.ddx.ddt.Ybar_p 
	X_q0.95 = q95.ddx.ddt.Ybar_p
	X_mean  = E.ddx.ddt.Ybar_p
	plot(t_,X_mean[,1],ylim=c(-3,3),cex=0,xlab="Time",ylab="Effects")
	for(j in 1:N)
	{
	  polygon(c(t_,rev(t_)),c(X_q0.05[,j],rev(X_q0.95[,j])),border=NA,col=rainbow(N,alpha=0.2)[j])
	  lines(t_,X_mean[,j],type="l",col=rainbow(N)[j])
	}
	lines(t_,rep(0,n_),lty=2)
	#
	## plot relative contributions
	X_q0.05 = q05.propC
	X_q0.95 = q95.propC
	X_mean  = E.propC
	plot(t_,X_mean[,1],ylim=c(0,1),cex=0,xlab="Time",ylab="Relative contribution")
	for(j in 1:N)
	{
	  polygon(c(t_,rev(t_)),c(X_q0.05[,j],rev(X_q0.95[,j])),border=NA,col=rainbow(N,alpha=0.2)[j])
	  lines(t_,X_mean[,j],type="l",col=rainbow(N)[j])
	}
	lines(t_,rep(0,n_),lty=2)
	#
	par(mfrow=c(1,1),mar=c(5,4,3,3)+0.1)
}

# ## plot mean effects and relative total contributions
# par(mfrow=c(2,1),mar=c(4,5,1,4),cex.lab=1)
# #
# .plot.diamond(meanE_ensemble, y_lab="Mean effect")
# .plot.diamond(totRelativeC_ensemble, y_lab="Total relative contribution")
# #
# par(mfrow=c(1,1))

# ## dynamical interaction plot 
# .plot.DIN(EMat,CMat,colnames(Y))

#
###
