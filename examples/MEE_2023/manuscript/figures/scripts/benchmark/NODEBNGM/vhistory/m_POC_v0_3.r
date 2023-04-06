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

## next:
##            - allow for different values of priors for each time series
##            - separate fitting script from visualisation script (?)
##            - extra visualisations (diamond plots, phase space, ...)

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
.heatMat <- function(YMat,lims=rbind(c(0,1),c(0,1)),labs=c("",""),main=c(""),axes=T)
{
  
  ## plot matrix
  #
  ## relative min and max to be matched to min and max color levels
  maxAbsMinMax <- max(abs(c(min(YMat),max(YMat))))
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
.plotMat <- function(lims,func,labs=c("",""),main=c(""),axes=T)
{
  YMat = .getMat(lims,func)
  .heatMat(YMat,lims,labs,main,axes)
}

#
###

#####################
## ARTIFICIAL DATA ##
#####################

## goal: simple LV with no interaction

## run ode
t_true     = seq(0,20,0.1)
Y_true_0   = c(1.0,1.0)
p          = c(1.0,0.5,1.0)
ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]*x[2],p[2]*x[1]*x[2]-p[3]*x[2]))
ddx.r_true = function(t,x,p) list(c(0,-p[2],p[2],0))
Ybar_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)

## create TS 
s = seq(0,20*10,5)
TS = Ybar_true[s,]
colnames(TS) = c("t","hare","lynx")
write.table(TS,file="data/TS_LV_1.csv",sep=";",row.names=FALSE)

## ground truth
Ybar_true     = Ybar_true[,2:3]
ddt.Ybar_true = t(apply(Ybar_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
ddx.r_true = t(apply(Ybar_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))

# ## visualise
# plot(model[,1],model[,2],col="red",type="l")
# lines(model[,1],model[,3],col="blue")


#
###

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
Ybar_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)

## create TS 
s = seq(1,20*10,5)
TS = Ybar_true[s,]
colnames(TS) = c("t","hare","lynx")
write.table(TS,file="data/TS_LV_2.csv",sep=";",row.names=FALSE)

## ground truth
Ybar_true     = Ybar_true[,2:3]
ddt.Ybar_true = t(apply(Ybar_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
ddx.r_true = t(apply(Ybar_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))

# ## visualise
# plot(TS[,1],TS[,2],col="red",type="l")
# lines(TS[,1],TS[,3],col="blue")

#
###

###################
## INITIATE NODE ##
###################

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_38.r")

## load data
TS = read.table("data/TS_LV_2.csv",sep=";",header=T)

## parameters
N_e    	 = 10      					# number of ensemble elements
W_o    	 = 10    					# number of nodes in observation model
W_p    	 = 10     					# number of nodes in process model
sd1_o  	 = 0.01   					# sd in likelihood of residuals of observation model
sd2_o  	 = c(rep(0.3 ,W_o),            
			 rep(0.01,W_o),
			 rep(0.2 ,W_o))	 		# sd in priors of observation model
sd1_p  	 = 0.01   	   				# sd in likelihood of residuals of process model
sd2_p  	 = 0.2  	   				# sd in priors of process model
alpha_i  = 2      	   				# factor of interpolation
alpha_o  = nrow(TS)-1  				# scaling factor for observation model

## initiate NODE
t     = TS[,1]                                    # vector of time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)    # normalise interpolated time steps
Y     = as.matrix(TS[,-1])                                   # matrix containing time series

#
###

pdf("results.pdf")

###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model

## standardise data
Y     = log(Y)
mean_ = mean(Y)
sd_   = sd(Y)
Y     = (Y-mean_)/sd_

## fit observation model
chainList = list()
for(i in 1:ncol(Y))
{
	message(paste("fitting: ",i,"/",ncol(Y),sep=""))
	Omega_o_chain = NULL
	for(k in 1:N_e)
	{
		Omega_o        = rnorm(W_o*3,0,sd2_o)
		logPost_o_0    = logPost_o(t,Y[,i],Omega_o,sd1_o,sd2_o)
		Omega_o        = argmax.logPost_o(t,Y[,i],Omega_o,sd1_o,sd2_o)
		logPost_o_f    = logPost_o(t,Y[,i],Omega_o,sd1_o,sd2_o)
		Omega_o_chain  = rbind(Omega_o_chain,c(logPost_o_f,Omega_o))

		message(paste("sample: ",k,"/",N_e," ",round(logPost_o_0,2)," --> ",round(logPost_o_f,2),sep=""))
	}
	chainList[[i]] = Omega_o_chain
	message("\n")
}

## keep best sample
Omega_o_MaP = list()
for(i in 1:ncol(Y))
{
	idx = which.max(chainList[[i]][,1])
	Omega_o_MaP[[i]] = chainList[[i]][idx,-1]
	message(max(chainList[[i]][,1]))
}
message("\n")

## evaluate observation model
Ybar_o_MaP     = NULL 
ddt.Ybar_o_MaP = NULL
for(i in 1:ncol(Y))
{
	Ybar_o_         = f_o.eval(t_,Omega_o_MaP[[i]])
	Ybar_o_MaP      = cbind(Ybar_o_MaP,Ybar_o_)
	ddt.Ybar_o_     = ddt.f_o.eval(t_,Omega_o_MaP[[i]])
	ddt.Ybar_o_MaP  = cbind(ddt.Ybar_o_MaP,ddt.Ybar_o_)
}

## evaluate ensemble
predictions_o = list()
for(i in 1:ncol(Y))
{
	Ybar_o     = t(apply(chainList[[i]][,-1],1,function(x) f_o.eval(t_,x)))
	Ybar_o     = exp(mean_ + sd_*Ybar_o)
	E.Ybar_o   = apply(Ybar_o,2,mean)
	q05.Ybar_o = apply(Ybar_o,2,quantile,p=0.05)
	q95.Ybar_o = apply(Ybar_o,2,quantile,p=0.95)

	ddt.Ybar_o     = t(apply(chainList[[i]][,-1],1,function(x) ddt.f_o.eval(t_,x)))
	ddt.Ybar_o     = sd_*Ybar_o*ddt.Ybar_o
	E.ddt.Ybar_o   = apply(ddt.Ybar_o,2,mean)
	q05.ddt.Ybar_o = apply(ddt.Ybar_o,2,quantile,p=0.05)
	q95.ddt.Ybar_o = apply(ddt.Ybar_o,2,quantile,p=0.95)

	predictions_o[[i]] = cbind(E.Ybar_o,q05.Ybar_o,q95.Ybar_o,E.ddt.Ybar_o,q05.ddt.Ybar_o,q95.ddt.Ybar_o)
}

## unscale data
Y              = exp(mean_ + sd_*Y)
Ybar_o_MaP     = exp(mean_ + sd_*Ybar_o_MaP)
ddt.Ybar_o_MaP = sd_*Ybar_o_MaP*ddt.Ybar_o_MaP

## visualise
par(mfrow=c(2,1))
for(i in 1:ncol(Y))
{
	##
	plot(t_true,Ybar_true[,i],lty=2,type="l")
	lines(t_,predictions_o[[i]][,"E.Ybar_o"],col="red")
	polygon(x=c(t_,rev(t_)),y=c(predictions_o[[i]][,"q05.Ybar_o"],rev(predictions_o[[i]][,"q95.Ybar_o"])),border=NA,col=adjustcolor("grey",alpha=0.5))
	lines(t_,Ybar_o_MaP[,i],type="l",col="red",lty=2)
	points(t,Y[,i])
	#
	##
	plot(t_true,ddt.Ybar_true[,i],lty=2,type="l")
	polygon(x=c(t_,rev(t_)),y=c(predictions_o[[i]][,"q05.ddt.Ybar_o"],rev(predictions_o[[i]][,"q95.ddt.Ybar_o"])),border=NA,col=adjustcolor("grey",alpha=0.5))
	lines(t_,predictions_o[[i]][,"E.ddt.Ybar_o"],type="l",col="red")
	lines(t_,ddt.Ybar_o_MaP[,i],type="l",col="red",lty=2)
}
par(mfrow=c(1,1))

## visualise periodic components
Ybar_o     = Ybar_o_MaP
ddt.Ybar_o = ddt.Ybar_o_MaP
for( i in 1:ncol(Y))
{
	Ybar_o_ = t(apply(matrix(Omega_o_MaP[[i]],ncol=3),1,function(x) f_o.eval(t_,x)))
	plot(t_,apply(Ybar_o_,2,sum),type="l")
	for(j in 1:W_o)
	{
		lines(t_,Ybar_o_[j,],col=rainbow(W_o)[j]) 
	}
}

## visualise prior and posterior distributions
par(mfrow=c(3,1))
for(i in 1:ncol(Y))
{
	for(j in 1:3)
	{
		tmp = as.vector(chainList[[i]][,-1][,1:W_o+W_o*(j-1)])
		sd2_o_ = sd2_o[1+W_o*(j-1)] 
		x = density(tmp)$x
		y = density(tmp)$y
		plot(x,1/sqrt(2*pi*sd2_o_^2)*exp(-(x^2)/(2*sd2_o_^2)),type="l",lty=2)
		lines(x,y)
	}
}
par(mfrow=c(1,1))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

# ##
# Ybar_o     = ode(y=Y_true_0,times=t_,func=ddt.Y_true,parms=p)[,2:3]
# ddt.Ybar_o = t(apply(Ybar_o,1,function(x)unlist(ddt.Y_true(NULL,x,p))))

##
Ybar_o     = Ybar_o_MaP
ddt.Ybar_o = ddt.Ybar_o_MaP

## standardise Ybar_o
mean_ = mean(Ybar_o)
sd_   = sd(Ybar_o)
std.Ybar_o = (Ybar_o-mean_)/sd_

## fit model
chainList = list()
for(i in 1:ncol(Y))
{
	message(paste("fitting: ",i,"/",ncol(Y),sep=""))
	Omega_p_chain = NULL
	for(k in 1:N_e)
	{
		Omega_p      = rnorm(W_p*(2+ncol(Y)),0,sd2_p)
		logPost_p_0  = logPost_p(std.Ybar_o,1/Ybar_o[,i]*ddt.Ybar_o[,i],Omega_p,sd1_p,sd2_p)
		Omega_p      = argmax.logPost_p(std.Ybar_o,1/Ybar_o[,i]*ddt.Ybar_o[,i],Omega_p,sd1_p,sd2_p)
		logPost_p_f  = logPost_p(std.Ybar_o,1/Ybar_o[,i]*ddt.Ybar_o[,i],Omega_p,sd1_p,sd2_p)
		Omega_p_chain  = rbind(Omega_p_chain,c(logPost_p_f,Omega_p))
		message(paste("sample: ",k,"/",N_e," ",round(logPost_p_0,2)," --> ",round(logPost_p_f,2),sep=""))
	}
	chainList[[i]] = Omega_p_chain
	message("\n")
}

## maximum a posteriori
Omega_p_MaP = list()
for(i in 1:ncol(Y))
{
	idx = which.max(chainList[[i]][,1])
	Omega_p_MaP[[i]] = chainList[[i]][idx,-1]
	message(max(chainList[[i]][,1]))
}
message("\n")

## evaluate best model
ddt.Ybar_p_MaP     = NULL
ddx.ddt.Ybar_p_MaP = NULL
for(i in 1:ncol(Y))
{
	ddt.Ybar_p_        = ddt.f_p.eval(std.Ybar_o,Omega_p_MaP[[i]])
	ddt.Ybar_p_MaP     = cbind(ddt.Ybar_p_MaP,ddt.Ybar_p_)
	ddx.ddt.Ybar_p_    = t(ddx.ddt.f_p.eval(std.Ybar_o,Omega_p_MaP[[i]]))
	ddx.ddt.Ybar_p_MaP = cbind(ddx.ddt.Ybar_p_MaP,ddx.ddt.Ybar_p_)
}

## evaluate ensemble
predictions_p = list()
for(i in 1:ncol(Y))
{
	ddt.Ybar_p     = t(apply(chainList[[i]][,-1],1,function(x) ddt.f_p.eval(std.Ybar_o,x)))
	E.ddt.Ybar_p   = apply(ddt.Ybar_p,2,mean)
	q05.ddt.Ybar_p = apply(ddt.Ybar_p,2,quantile,p=0.05)
	q95.ddt.Ybar_p = apply(ddt.Ybar_p,2,quantile,p=0.95)

	ddx.ddt.Ybar_p     = t(apply(chainList[[i]][,-1],1,function(x) ddx.ddt.f_p.eval(std.Ybar_o,x)))
	ddx.ddt.Ybar_p     = 1/sd_*ddx.ddt.Ybar_p
	E.ddx.ddt.Ybar_p   = t(matrix(apply(ddx.ddt.Ybar_p,2,mean),ncol=length(t_)))
	q05.ddx.ddt.Ybar_p = t(matrix(apply(ddx.ddt.Ybar_p,2,quantile,p=0.05),ncol=length(t_)))
	q95.ddx.ddt.Ybar_p = t(matrix(apply(ddx.ddt.Ybar_p,2,quantile,p=0.95),ncol=length(t_)))

	predictions_p[[i]] = list(E.ddt.Ybar_p,q05.ddt.Ybar_p,q95.ddt.Ybar_p,E.ddx.ddt.Ybar_p,q05.ddx.ddt.Ybar_p,q95.ddx.ddt.Ybar_p)
}

## unscale data
ddx.ddt.Ybar_p_MaP = 1/sd_*ddx.ddt.Ybar_p_MaP

## visualise predictions
par(mfrow=c(2,1))
for(i in 1:ncol(Y))
{
	##
	plot(t_true,1/Ybar_true[,i]*ddt.Ybar_true[,i],lty=2,type="l")
	polygon(x=c(t_,rev(t_)),y=c(predictions_p[[i]][[2]],rev(predictions_p[[i]][[3]])),border=NA,col=adjustcolor("grey",alpha=0.5))
	lines(t_,1/Ybar_o[,i]*ddt.Ybar_o[,i],type="l",col="black")
	lines(t_,ddt.Ybar_p_MaP[,i],type="l",col="red")
	lines(t_,predictions_p[[i]][[1]],type="l",col="red",lty=2)

	##
	plot(t_,ddx.ddt.Ybar_p_MaP[,1+(i-1)*2],ylim=c(-2,2),col="white")
	for(j in 1:ncol(Y))
	{
		polygon(x=c(t_,rev(t_)),y=c(predictions_p[[i]][[5]][,j],rev(predictions_p[[i]][[6]][,j])),border=NA,col=adjustcolor("grey",alpha=0.5))
		lines(t_,predictions_p[[i]][[4]][,j],type="l",col=j)
		lines(t_true,ddx.r_true[,j+(i-1)*2],type="l",col=j,lty=2)
	}	
}
par(mfrow=c(1,1))

## visualise NODE approx
par(mfrow=c(2,2))
for(i in 1:ncol(Y))
{
	func = function(x)
	{	
		x = x*sd_ + mean_
		r = 1/x[i]*unlist(ddt.Y_true(NULL,x,p))[i]
		return(r)
	}
	.plotMat(rbind(c(-3,3),c(-3,3)),func)
	y = (std.Ybar_o - min(std.Ybar_o))/(max(std.Ybar_o)-min(std.Ybar_o))
	#
	func = function(x) ddt.f_p(x,Omega_p_MaP[[i]])
	.plotMat(rbind(c(-3,3),c(-3,3)),func)
}

#
###

dev.off()
