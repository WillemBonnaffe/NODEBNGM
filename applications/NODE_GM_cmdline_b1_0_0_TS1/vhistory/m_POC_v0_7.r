###########
## m_POC ##
###########

## goal: perform assessment of quality of inference of NODE by gradient matching 

## update:
## 03-03-2022 - created v0_0 from m_NODE_GM_v0_17.r
## 07-03-2022 - created v0_1
##            - re-introduced scaling of variables to help with the fitting 
##            - created v0_2
##            - tested DEMCpp and DEMCO
## 08-03-2022 - created v0_3
##            - coded anchor sampling
##            - added specific priors
##            - re-introduced FFT as benchmark
## 10-03-2022 - created v0_5
##            - cleaned code (removed unused sections)
## 11-03-2022 - created v0_6
##            - implemented FFT and DFT
##            - added tri-trophic time series 
## 14-03-2022 - created v0_7
##            - compressed and cleaned code

## next:
##            - individual scaling descaling

##############
## INITIATE ##
##############

## 
library(deSolve)

#
###

#######################
## UTILITY FUNCTIONS ##
#######################

## goal: functions for visualisation (figure 3 in the paper)

## .getmat ##
## goal: return matrix of function evaluated at grid points in variable limits 
## args: 
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
.getMat <- function(lims,func)
{
  
  ## compute 2D matrix
  res <- 100
  X <- apply(lims,1,function(X){seq(X[1],X[2],(X[2]-X[1])/res)}) # grid of pairwise combinations of x and y
  #
  ## fill output matrix - function value at each c(x,y)
  YMat <- matrix(0,ncol=res+1,nrow=res+1)
  for(i in 1:length(X[,1]))
  {
    for(j in 1:length(X[,2]))
    {
      YMat[i,j] <- func(c(X[i,1],X[j,2]))
    }
  }
  
  ## return matrix
  return(YMat)

}

## .heatmat ##
## goal: heatmap function for plotting 2D matrix
## args: 
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
  levels <- seq(-maxAbsMinMax,maxAbsMinMax,2*maxAbsMinMax/1000)
  colorLevels <- rev(rainbow(1000,start=0,end=1,alpha=0.5))
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
## args: 
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

## DFT
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

#
###

pdf("results.pdf")

# #####################
# ## ARTIFICIAL DATA ##
# #####################
# 
# ## goal: simple LV with no interaction
# 
# ## run ode
# t_true     = seq(0,20,0.1)
# Y_true_0   = c(1.0,1.0)
# p          = c(1.0,0.5,1.0)
# ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]*x[2],p[2]*x[1]*x[2]-p[3]*x[2]))
# ddx.r_true = function(t,x,p) list(c(0,-p[2],p[2],0))
# Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)
# 
# ## create TS 
# s = seq(0,20*10,5)
# TS = Yhat_true[s,]
# colnames(TS) = c("t","hare","lynx")
# write.table(TS,file="data/TS_LV_1.csv",sep=";",row.names=FALSE)
# 
# ## ground truth
# Yhat_true     = Yhat_true[,2:3]
# ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
# ddx.r_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))
# 
# #
# ###
 
#####################
## ARTIFICIAL DATA ##
#####################

## goal: simple LV with no interaction

## run ode
t_true     = seq(0,20,0.1)
Y_true_0   = c(1.0,1.0)
p          = c(1.0,1.5,0.25,1.0)
ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]/(1+p[2]*p[3]*x[1])*x[2],p[2]*x[1]/(1+p[2]*p[3]*x[1])*x[2]-p[4]*x[2]))
ddx.r_true = function(t,x,p) list(c(+p[2]*x[2]*p[2]*p[3]/(1+p[2]*p[3]*x[1])^2,
								    -p[2]/(1+p[2]*p[3]*x[1]),
                                     p[2]/(1+p[2]*p[3]*x[1])^2,
                                     0))
Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)

## create TS
s = seq(1,20*10,5)
TS = Yhat_true[s,]
colnames(TS) = c("t","hare","lynx")
write.table(TS,file="data/TS_LV_2.csv",sep=";",row.names=FALSE)

## ground truth
Yhat_true     = Yhat_true[,2:3]
ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
ddx.r_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))

## visualise
plot(TS[,1],TS[,2],col="red",type="l")
lines(TS[,1],TS[,3],col="blue")

#
###

# #####################
# ## ARTIFICIAL DATA ##
# #####################
# 
# ## goal: simple LV with no interaction
# 
# ## run ode
# t_true     = seq(0,20,0.1)
# Y_true_0   = c(1.0,1.0,0.1)
# p          = c(1.0,0.5,0.25,0.4,1.0,1.0)
# ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]*x[2]-p[3]*x[1]*x[3],
# 								              p[2]*x[1]*x[2]-p[4]*x[2]*x[3]-p[5]*x[2],
# 									          p[3]*x[1]*x[3]+p[4]*x[2]*x[3]-p[6]*x[3]))
# ddx.r_true = function(t,x,p) list(c(0,-p[2],-p[3],p[2],0,-p[4],p[3],p[4],0))
# Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)
# 
# ## create TS 
# s = seq(0,20*10,5)
# TS = Yhat_true[s,]
# colnames(TS) = c("t","R","G","B")
# write.table(TS,file="data/TS_LV_3.csv",sep=";",row.names=FALSE)
# 
# ## ground truth
# Yhat_true     = Yhat_true[,2:4]
# ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
# ddx.r_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))
# 
# ## visualise
# plot(TS[,1],TS[,2],col="red",type="l",ylim=c(0,max(TS[,-1])))
# lines(TS[,1],TS[,3],col="green")
# lines(TS[,1],TS[,4],col="blue")
# 
# #
# ###

###################
## INITIATE NODE ##
###################

## goal: initiate the NODE 

## load data
# TS = read.table("data/TS_LV_3.csv",sep=";",header=T)

## initiate TS 
t     = 1:nrow(TS) # TS[,1]                                    # vector of time steps
Y     = as.matrix(TS[,-1])                                   # matrix containing time series
alpha_i  = 2      	   				                          # factor of interpolation

#
###
 
###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model with DFT and FFT

##
dt    = diff(TS[1:2,1])
t     = 1:nrow(TS) # TS[,1]                       # vector of time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)    # normalise interpolated time steps

## standardise data
Y     = log(Y)
mean_ = mean(Y)
sd_   = sd(Y)
Y     = (Y-mean_)/sd_

## initiate
up = 2
n  = 10
f  = Y[,1]
t  = seq(1,length(f))
nt = seq(1,length(f),1/up)

## compute DFT
Yhat_o_DFT     = NULL
ddt.Yhat_o_DFT = NULL
for(i in 1:ncol(Y))
{
	Yhat_o_DFT = cbind(Yhat_o_DFT,DFT(Y[,i],t,nt,n))
	ddt.Yhat_o_DFT = cbind(ddt.Yhat_o_DFT,dDFT(Y[,i],t,nt,n))
}

## unscale data
t  = (t-min(t))/(max(t)-min(t))*(max(TS[,1])-min(TS[,1]))+min(TS[,1])
nt = (nt-min(nt))/(max(nt)-min(nt))*(max(TS[,1])-min(TS[,1]))+min(TS[,1])
t_ = (t_-min(t_))/(max(t_)-min(t_))*(max(TS[,1])-min(TS[,1]))+min(TS[,1])

## unscale data
Y              = exp(mean_ + sd_*Y)
Yhat_o_DFT     = exp(mean_ + sd_*Yhat_o_DFT)
ddt.Yhat_o_DFT = 1/dt*sd_*Yhat_o_DFT*ddt.Yhat_o_DFT
#
# Y              = mean_ + sd_*Y
# Yhat_o_DFT     = mean_ + sd_*Yhat_o_DFT
# ddt.Yhat_o_DFT = 1/dt*sd_*ddt.Yhat_o_DFT

## visualise
par(mfrow=c(2,1))
#
for(i in 1:ncol(Y))
{
	plot(t,TS[,i+1])
	lines(nt,Yhat_o_DFT[,i],col="blue",lty=1)
	#
	plot(t_true,ddt.Yhat_true[,i],type="l",lty=2,ylim=c(-1,1)*3)
	lines(nt,ddt.Yhat_o_DFT[,i],col="blue",lty=1)
}
#
par(mfrow=c(1,1))

#
###

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load data
# TS = read.table("data/TS_LV_2.csv",sep=";",header=T)

## load NODE functions
source("f_NODE_GM_v0_38.r")

#
###

###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model

## parameters
N_e    	= 10    					                            # number of ensemble elements
W_o    	= 10    					                            # number of nodes in observation model
sd1_o  	= 0.01   					                            # sd in likelihood of residuals of observation model
sd2_o  	= c(rep(0.2,W_o),rep(0.1,W_o),rep(0.1,W_o))	 		# sd in priors of observation model
alpha_i = 2      	   											# factor of interpolation
results_o = list()
for (i in 1:(ncol(TS)-1))
{
	## 
	message(paste("fitting: ",i,"/",ncol(TS)-1,sep=""))

	## initiate time series
	dt    = diff(TS[1:2,1])
	t     = 1:nrow(TS) 		 	                      # vector of time steps
	t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)    # vector of interpolated time steps
	Y     = as.matrix(TS[,i+1])                       # matrix containing time series
	
	## standardise/transform data
	Y     = log(Y)                                    # log facilitates fitting
	mean_ = mean(Y)
	sd_   = sd(Y)
	Y     = (Y-mean_)/sd_
	
	## fit observation model
	Omega_o_chain = NULL
	for(k in 1:N_e)
	{
		Omega_o_0      = rnorm(W_o*3,0,sd2_o)
		Omega_o_f      = argmax.logPost_o(t,Y,Omega_o_0,sd1_o,sd2_o)
		logPost_o_0    = logPost_o(t,Y,Omega_o_0,sd1_o,sd2_o)
		logPost_o_f    = logPost_o(t,Y,Omega_o_f,sd1_o,sd2_o)
		Omega_o_chain  = rbind(Omega_o_chain,c(logPost_o_f,Omega_o_f))
		message(paste(k,"/",N_e,"\t",
			round(logPost_o_0,2),"\t","-->","\t",
			round(logPost_o_f,2),sep=""))
	}
	message("")
	
	## maximum a posteriori (MaP)
	idx            = which.max(Omega_o_chain[,1])
	Omega_o_MaP    = Omega_o_chain[idx,-1]
	MaP.Yhat_o     = f_o.eval(t_,Omega_o_MaP)
	MaP.ddt.Yhat_o = ddt.f_o.eval(t_,Omega_o_MaP)
	message(paste("MaP: ",max(Omega_o_chain[idx,1],sep="")))
	message("\n")
	
	## evaluate ensemble
	Yhat_o         = t(apply(Omega_o_chain[,-1],1,function(x) f_o.eval(t_,x)))
	E.Yhat_o       = apply(Yhat_o,2,mean)
	q05.Yhat_o     = apply(Yhat_o,2,quantile,p=0.05)
	q95.Yhat_o     = apply(Yhat_o,2,quantile,p=0.95)
	ddt.Yhat_o     = t(apply(Omega_o_chain[,-1],1,function(x) ddt.f_o.eval(t_,x)))
	E.ddt.Yhat_o   = apply(ddt.Yhat_o,2,mean)
	q05.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.05)
	q95.ddt.Yhat_o = apply(ddt.Yhat_o,2,quantile,p=0.95)
	
	## de-scale time [1,N] -> [0,max(t)]
	t  = (t-min(t))/(max(t)-min(t))*(max(TS[,1])-min(TS[,1]))+min(TS[,1])
	t_ = (t_-min(t_))/(max(t_)-min(t_))*(max(TS[,1])-min(TS[,1]))+min(TS[,1])
	
	## de-scale data
	Y              = exp(mean_ + sd_*Y)
	MaP.Yhat_o     = exp(mean_ + sd_*MaP.Yhat_o)
	E.Yhat_o       = exp(mean_ + sd_*E.Yhat_o)
	q05.Yhat_o     = exp(mean_ + sd_*q05.Yhat_o)
	q95.Yhat_o     = exp(mean_ + sd_*q95.Yhat_o)
	MaP.ddt.Yhat_o = 1/dt*sd_*MaP.Yhat_o*MaP.ddt.Yhat_o
	E.ddt.Yhat_o   = 1/dt*sd_*MaP.Yhat_o*E.ddt.Yhat_o
	q05.ddt.Yhat_o = 1/dt*sd_*MaP.Yhat_o*q05.ddt.Yhat_o
	q95.ddt.Yhat_o = 1/dt*sd_*MaP.Yhat_o*q95.ddt.Yhat_o
	
	## store
	results_o[[i]] = list(t_,MaP.Yhat_o,MaP.ddt.Yhat_o)
	names(results_o[[i]]) = c("t_","MaP.Yhat_o","MaP.ddt.Yhat_o")

	## visualise
	par(mfrow=c(2,1))
	#
	## plot density
	plot(t_,rep(0,length(t_)),ylim=c(0,1)*max(Y),type="l",lty=3)
	points(t,Y)
	lines(t_true,Yhat_true[,i],lty=2)
	polygon(x=c(t_,rev(t_)),y=c(q05.Yhat_o,rev(q95.Yhat_o)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.Yhat_o,col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_,MaP.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i],lty=2)
	# lines(nt,Yhat_o_DFT[,i],type="l",col="blue",lty=2)
	#
	## plot dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*4,type="l",lty=3)
	lines(t_true,ddt.Yhat_true[,i],lty=2)
	polygon(x=c(t_,rev(t_)),y=c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.ddt.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_,MaP.ddt.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i],lty=2)
	# lines(nt,ddt.Yhat_o_DFT[,i],type="l",col="blue",lty=2)
	#
	par(mfrow=c(1,1))
	
	## visualise periodic components
	Yhat_o_      = t(apply(matrix(Omega_o_MaP,ncol=3),1,function(x) f_o.eval(t_,x)))
	sum.Yhat_o_  = apply(Yhat_o_,2,sum)
	plot(t_,sum.Yhat_o_,ylim=c(min(Yhat_o_),max(Yhat_o_)),cex=0)
	for(j in 1:W_o)
	{
		lines(t_,Yhat_o_[j,],col=rainbow(W_o)[j])
	}
	
	## visualise prior and posterior distributions
	par(mfrow=c(3,1))
	for(j in 1:3)
	{
		chain_ = as.vector(Omega_o_chain[,-1][,1:W_o+W_o*(j-1)])
		sd2_o_ = sd2_o[1+W_o*(j-1)]
		x      = density(chain_)$x
		y      = density(chain_)$y
		pd     = 1/sqrt(2*pi*sd2_o_^2)*exp(-(x^2)/(2*sd2_o_^2))
		plot(x,pd,type="l",lty=2,ylim=c(0,max(c(y,pd))))
		lines(x,y)
	}
	par(mfrow=c(1,1))

}

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## parameters
N_e    	 = 10     					# number of ensemble elements
W_p    	 = 10     					# number of nodes in process model
sd1_p  	 = 0.01   	   				# sd in likelihood of residuals of process model
sd2_p  	 = 0.1  	   				# sd in priors of process model

## prepare data
Yhat_o     = NULL 
ddt.Yhat_o = NULL
for(i in 1:(ncol(TS)-1))
{
	Yhat_o     = cbind(Yhat_o,    results_o[[i]]$MaP.Yhat_o)
	ddt.Yhat_o = cbind(ddt.Yhat_o,results_o[[i]]$MaP.ddt.Yhat_o)
}
#
## or use DFT
# Yhat_o     = Yhat_o_DFT 
# ddt.Yhat_o = ddt.Yhat_o_DFT 
#
## or ground truth
# Yhat_o     = ode(y=Y_true_0,times=t_,func=ddt.Y_true,parms=p)[,2:3]
# ddt.Yhat_o = t(apply(Yhat_o,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
#
Y = 1/Yhat_o*ddt.Yhat_o

## standardise predictive variables
mean_x = apply(Yhat_o,2,mean)
sd_x   = apply(Yhat_o,2,sd)
Yhat_o = t((t(Yhat_o)-mean_x)/sd_x)

i = 1
results_p = list()
for(i in 1:(ncol(TS)-1))
{
	## 
	message(paste("fitting: ",i,"/",ncol(TS)-1,sep=""))

	## standardise response
	r_o    = Y[,i]
	mean_y = mean(r_o)
	sd_y   = sd(r_o) 
	r_o    = (r_o-mean_y)/sd_y

	## fit model
	Omega_p_chain = NULL
	for(k in 1:N_e)
	{
		Omega_p_0      = rnorm(W_p*(2+ncol(TS)-1),0,sd2_p)
		Omega_p_f      = argmax.logPost_p(Yhat_o,r_o,Omega_p_0,sd1_p,sd2_p)
		logPost_p_0    = logPost_p(Yhat_o,r_o,Omega_p_0,sd1_p,sd2_p)
		logPost_p_f    = logPost_p(Yhat_o,r_o,Omega_p_f,sd1_p,sd2_p)
		Omega_p_chain  = rbind(Omega_p_chain,c(logPost_p_f,Omega_p_f))
		message(paste(k,"/",N_e,"\t",
				round(logPost_p_0,2),"\t","-->","\t",
				round(logPost_p_f,2),sep=""))
	}
	message("")
	
	## maximum a posteriori
	idx         = which.max(Omega_p_chain[,1])
	Omega_p_MaP = Omega_p_chain[idx,-1]
	message(Omega_p_chain[idx,1])
	message("\n")
	
	## evaluate best model
	MaP.r_p     = ddt.f_p.eval(Yhat_o,Omega_p_MaP)
	MaP.ddx.r_p = t(ddx.ddt.f_p.eval(Yhat_o,Omega_p_MaP))
	
	## evaluate ensemble
	r_p         = t(apply(Omega_p_chain[,-1],1,function(x) ddt.f_p.eval(Yhat_o,x)))
	E.r_p       = apply(r_p,2,mean)
	q05.r_p     = apply(r_p,2,quantile,p=0.05)
	q95.r_p     = apply(r_p,2,quantile,p=0.95)
	ddx.r_p     = t(apply(Omega_p_chain[,-1],1,function(x) ddx.ddt.f_p.eval(Yhat_o,x)))
	E.ddx.r_p   = t(matrix(apply(ddx.r_p,2,mean),ncol=length(t_)))
	q05.ddx.r_p = t(matrix(apply(ddx.r_p,2,quantile,p=0.05),ncol=length(t_)))
	q95.ddx.r_p = t(matrix(apply(ddx.r_p,2,quantile,p=0.95),ncol=length(t_)))
	Geber_p     = t(as.vector(t(Yhat_o))*t(ddx.r_p))
	E.Geber_p   = t(matrix(apply(Geber_p,2,mean),ncol=length(t_)))
	q05.Geber_p = t(matrix(apply(Geber_p,2,quantile,p=0.05),ncol=length(t_)))
	q95.Geber_p = t(matrix(apply(Geber_p,2,quantile,p=0.95),ncol=length(t_)))
	
	## unscale data
	r_o         = mean_y + sd_y * r_o 
	MaP.r_p     = mean_y + sd_y * MaP.r_p
	E.r_p       = mean_y + sd_y * E.r_p
	q05.r_p     = mean_y + sd_y * q05.r_p
	q95.r_p     = mean_y + sd_y * q95.r_p
	MaP.ddx.r_p = sd_y * t(1/sd_x * t(MaP.ddx.r_p))
	E.ddx.r_p   = sd_y * t(1/sd_x * t(E.ddx.r_p))
	q05.ddx.r_p = sd_y * t(1/sd_x * t(q05.ddx.r_p))
	q95.ddx.r_p = sd_y * t(1/sd_x * t(q95.ddx.r_p))
	
	## store
	results_p[[i]] = list(E.r_p,q05.r_p,q95.r_p,E.ddx.r_p,q05.ddx.r_p,q95.ddx.r_p)
	names(results_p[[i]]) = c("E.r_p","q05.r_p","q95.r_p","E.ddx.r_p","q05.ddx.r_p","q95.ddx.r_p")
	
	## visualise predictions
	par(mfrow=c(3,1))
	#
	## plot y
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*2,type="l",lty=3)
	points(t_,r_o)
	polygon(x=c(t_,rev(t_)),y=c(q05.r_p,rev(q95.r_p)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.r_p,col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_true,1/Yhat_true[,i]*ddt.Yhat_true[,i],lty=2)
	#
	## plot derivative
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*1,type="l",lty=3)
	for(j in 1:(ncol(TS)-1))
	{
		polygon(x=c(t_,rev(t_)),y=c(q05.ddx.r_p[,j],rev(q95.ddx.r_p[,j])),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[j])
		lines(t_,E.ddx.r_p[,j],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j])
		lines(t_true,ddx.r_true[,j+(i-1)*(ncol(TS)-1)],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j],lty=2)
	}
	#
	## plot Geber
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*2,type="l",lty=3)
	for(j in 1:(ncol(TS)-1))
	{
		polygon(x=c(t_,rev(t_)),y=c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[j])
		lines(t_,E.Geber_p[,j],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j])
	}
	#
	par(mfrow=c(1,1))
	
	## visualise prior and posterior distributions
	par(mfrow=c(3,1))
	for(j in 1:3)
	{
		tmp = as.vector(Omega_p_chain[,-1][,1:W_p+W_p*(j-1)])
		# sd2_p_ = sd2_[1+W_p[i]*(j-1)]
		x = density(tmp)$x
		y = density(tmp)$y
		pd = 1/sqrt(2*pi*sd2_p^2)*exp(-(x^2)/(2*sd2_p^2))
		plot(x,pd,type="l",lty=2,ylim=c(0,max(c(y,pd))))
		lines(x,y)
	}
	par(mfrow=c(1,1))
	
	# ## visualise NODE approx
	# par(mfrow=c(2,2))
	# #
	# func = function(x)
	# {
	# 	x = mean_x  + sd_x * x 
	# 	r = 1/x[i]*unlist(ddt.Y_true(NULL,x,p))[i]
	# 	return(r)
	# }
	# .plotMat(rbind(c(-3,3),c(-3,3)),func,maxAbsMinMax=5)
	# #
	# func = function(x) mean_y + sd_y * ddt.f_p(x,Omega_p_MaP)
	# .plotMat(rbind(c(-3,3),c(-3,3)),func,maxAbsMinMax=5)
	# #
	# par(mfrow=c(1,1))
}

#
###

dev.off()
