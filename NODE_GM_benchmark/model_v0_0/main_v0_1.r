############
## main.r ##
############

## goal: main template to fit a 3D NODE model to a time series

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update
## 12-12-2020 - created template for NODER v1.0.0
## 06-06-2022 - created v0_1
##            - adapted template to fit time series with 3 variables

## to do next:
## 09-06-2022 - implement de-standardisation of variables and Jacobian after fitting

##############
## INITIATE ##
##############

## import NODE functions
source("https://raw.githubusercontent.com/WillemBonnaffe/NODER/Branch_1/versions/v1.0/NODER_v1.0.1.r")

#
###

###################
## NODE INITIATE ##
###################

## load time series
TS = read.table("data/TS_3DLV.csv",sep=";",header=T)

## format time series
colnames(TS) = c("t","R","G","B") # compact names
TS = TS[20:50,]
# TS_ = TS.format(TS,c(2,3,4))    # log transform and standardise
TS[,1] = TS[,1] - min(TS[,1])
TS_  = TS

## initiate model
M0 = NODE.init(TS=TS_,
			   SLPinputStateList=list(c(2,3,4),c(2,3,4),c(2,3,4)),
			   SLPftypeVect=c(5,5,5)
			   )

## or load previous model
# load("model.RData")

## user defined parameters for the normal Bayesian model
sd_lik = c(0.1,0.1,0.1)
sd_prior_network = c(10.0,10.0,10.0)
#
## initiate sd_lik
sd_lik = rep(sd_lik,rep(M0$n,M0$N-1))
#
## initiate mu prior
mu_prior_states  = as.vector(t(M0$TS[1,-1])) # set mu_prior_states to first observation
mu_prior_network = rep(0,M0$D-(M0$N-1))
mu_prior         = c(mu_prior_states,mu_prior_network)
#
## initiate sd prior
sd_prior_states  = as.vector(t(apply(M0$TS[,-1],2,sd))) # set sd_prior_states to sd of time series
sd_prior_network = rep(sd_prior_network,M0$dimVect[-1]) 
sd_prior         = c(sd_prior_states,sd_prior_network)

## custom NODE Ydot
M0$D = M0$D*2
NODE.Ydot = function(state,OmegaList)
{	
  Ydot    = 1 # i = 1 is time => dY/dt[1] = 1
  Ydot[2] = (SLP(state[c(2,3,4)],OmegaList[[2]][1:(M0$D/2)],c(3,10,1),f_sigmaList[[1]]) + SLP(state[c(2,3,4)],OmegaList[[2]][(M0$D/2+1):M0$D],c(3,10,1),f_sigmaList[[5]]))*state[2]
  Ydot[3] = (SLP(state[c(2,3,4)],OmegaList[[3]][1:(M0$D/2)],c(3,10,1),f_sigmaList[[1]]) + SLP(state[c(2,3,4)],OmegaList[[3]][(M0$D/2+1):M0$D],c(3,10,1),f_sigmaList[[5]]))*state[3]
  Ydot[4] = (SLP(state[c(2,3,4)],OmegaList[[4]][1:(M0$D/2)],c(3,10,1),f_sigmaList[[1]]) + SLP(state[c(2,3,4)],OmegaList[[4]][(M0$D/2+1):M0$D],c(3,10,1),f_sigmaList[[5]]))*state[4]
  return(Ydot)
}
M0$NODE.Ydot = NODE.Ydot

## custom NODE Jacobian
NODE.J = function(state,OmegaList)
{	
  J = matrix(0,N,N) # i = 1 is time => J[1,] = d(dY/dt_1)/dY_1, d(dY/dt_1)/dY[2], ... = 0
  J[2,c(2,3,4)] = dSLP(state[c(2,3,4)],OmegaList[[2]][1:(M0$D/2)],c(3,10,1),df_sigmaList[[1]]) + dSLP(state[c(2,3,4)],OmegaList[[2]][(M0$D/2+1):M0$D],c(3,10,1),df_sigmaList[[5]])
  J[3,c(2,3,4)] = dSLP(state[c(2,3,4)],OmegaList[[3]][1:(M0$D/2)],c(3,10,1),df_sigmaList[[1]]) + dSLP(state[c(2,3,4)],OmegaList[[3]][(M0$D/2+1):M0$D],c(3,10,1),df_sigmaList[[5]])
  J[4,c(2,3,4)] = dSLP(state[c(2,3,4)],OmegaList[[4]][1:(M0$D/2)],c(3,10,1),df_sigmaList[[1]]) + dSLP(state[c(2,3,4)],OmegaList[[4]][(M0$D/2+1):M0$D],c(3,10,1),df_sigmaList[[5]])
  return(J)
}

#
###

##############
## NODE FIT ##
##############

# ## prefit with LASSO/Ridge regression
# # of dYdt operator on raw time differences from time series
# error.prefit = function(Omega)
# {
#   .error.LR(Y      =  TS.timedifferentiate(M0$TS)[,-1],
#             X      =  M0$TS[-M0$n,],
#             f      =  function(X,param){t(apply(X,1,function(x)M0$.Ydot(x,param)[-1]))},
#             param  = Omega,
#             lambda = 0.001,
#             alpha  = 1
#             )[[1]]
# }
# Omega = rnorm(M0$D-(M0$N-1),0,0.1); print(error.prefit(Omega))
# optList = .fit(Omega,error.prefit,20)
# Omega = optList$par

## fit
error.fit = function(Omega)
{
  -.error.normalBayesian(Y         = M0$TS[,-1],
                         X         = M0$TS[,1],
                         f         = function(X,param){M0$.Ybar.ode(X,param[M0$i_init],param[-M0$i_init])[,-1]},
                         param     = Omega,
                         sd_lik    = sd_lik,
                         mu_prior  = mu_prior,
                         sd_prior  = sd_prior
                         )[[1]]
} 
M0$Theta = c(mu_prior_states,rnorm(M0$D-(M0$N-1),0,0.01)) # Theta = rnorm(M0$D,0,0.01)
# load("model.RData")
error.fit(M0$Theta)
optList = .fit(M0$Theta,error.fit,1000)
M0$Theta = optList$par

## estimate error with ABC
chainList = .ABC(M0$Theta,error.fit,100,noise=0.001,threshold=0.5)
M0$chain = chainList$chain

## save model
save(M0,file="model.RData")

#
###

###################
## NODE ANALYSIS ##
###################

pdf("results.pdf")

## visualise dynamics
t = seq(min(M0$TS[,1]),max(M0$TS[,1]),0.1)
Y_0 = M0$Theta[M0$i_init]
Omega = M0$Theta[-M0$i_init]
Ybar = M0$.Ybar.ode(t,Y_0,Omega)
.plot.fit(M0$TS,Ybar)

## check residuals
res = M0$TS[,-1] - M0$.Ybar.ode(M0$TS[,1],Y_0,Omega)[,-1]
.plot.residualCheck1(res)
.plot.residualCheck2(res)
.plot.residualCheck3(res)
r2 = 1-sd(t(res))^2/sd(t(M0$TS[,-1]))^2

## compute approximated posterior distribution
t = seq(0,30,1)
#
## states
Ybar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Ybar.ode(t,x[M0$i_init],x[-M0$i_init]))
Ybar.q05 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.05))),ncol=M0$N)
Ybar.q50 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.50))),ncol=M0$N)
Ybar.q95 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.95))),ncol=M0$N)
#
## jacobian
Jbar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Jbar.ode(t,x[M0$i_init],x[-M0$i_init]))
Jbar.q05 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.05)))),nrow=length(t))
Jbar.q50 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.50)))),nrow=length(t))
Jbar.q95 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.95)))),nrow=length(t))
#
## contributions
Cbar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Cbar.ode(t,x[M0$i_init],x[-M0$i_init]))
Cbar.q05 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.05)))),nrow=length(t))
Cbar.q50 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.50)))),nrow=length(t))
Cbar.q95 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.95)))),nrow=length(t))

## plot fit, effects, contributions
layout(mat=rbind(c(1,4),c(2,5),c(3,6)))
colVect = rainbow(M0$N)
#
for(i in 2:M0$N)
{
	## states
	X.q05 = Ybar.q05
	X.q50 = Ybar.q50
	X.q95 = Ybar.q95
	plot(M0$TS[,1],M0$TS[,i],col=colVect[i],xlim=c(min(t),max(t)),ylim=c(-1,1)*3,xlab="t",ylab=colnames(M0$TS)[i])
	polygon(x=c(t,rev(t)),y=c(X.q05[,i],rev(X.q95[,i])),border=NA,col=adjustcolor(colVect[i],alpha=0.2))
	lines(t,X.q50[,i],col=colVect[i])
	#
	## jacobian
	X.q05 = Jbar.q05
	X.q50 = Jbar.q50
	X.q95 = Jbar.q95
	plot(x=c(min(t),max(t)),y=c(0,0),lty=2,ylim=c(-1,1)*1.0,cex=0,xlab="t",ylab="Effects")
	for(j in 1:M0$N)
	{
		k = (j-1)*(M0$N) + i
		polygon(x=c(t,rev(t)),y=c(X.q05[,k],rev(X.q95[,k])),border=NA,col=adjustcolor(colVect[j],alpha=0.2))
		lines(t,X.q50[,k],col=colVect[j])
	}
	#
	## contributions
	X.q05 = Cbar.q05
	X.q50 = Cbar.q50
	X.q95 = Cbar.q95
	plot(x=c(min(t),max(t)),y=c(0,0),lty=2,ylim=c(-1,1)*1.0,cex=0,xlab="t",ylab="Contributions")
	for(j in 1:M0$N)
	{
		k = (j-1)*(M0$N) + i
		polygon(x=c(t,rev(t)),y=c(X.q05[,k],rev(X.q95[,k])),border=NA,col=adjustcolor(colVect[j],alpha=0.2))
		lines(t,X.q50[,k],col=colVect[j])
	}
}
#
par(mfrow=c(1,1))

## average Jacobian and contributions
Jbar.q50.mean = matrix(apply(Jbar.q50,2,mean),ncol=M0$N)[-1,-1]
Cbar.q50.mean = matrix(apply(abs(Cbar.q50),2,mean),ncol=M0$N)[-1,-1]
#
## dynamical interaction network (DIN) plot
.plot.DIN(Jbar.q50.mean,Cbar.q50.mean,labels=colnames(M0$TS)[-1])

# ## visualise dYdt functions
# i=2
# for(i in c(2,3))
# {
# 	dYdt = function(x){M0$.Ydot(c(0,x[1],x[2]),M0$Theta[-M0$i_init])[i]}
# 	.plotMat(lims=rbind(c(-2,2),c(-2,2)),func=dYdt)
# }
# 
# ## pairwise phase space (PPS) plot
# dYdt = function(x){M0$.Ydot(c(0,x[1],x[2]),M0$Theta[-M0$i_init])}
# .plot.PPS(2,3,c(2,3),dYdt,rbind(c(-2,2),c(-2,2)))
# 
# ## find equilibrium
# Y_eq = rnorm(M0$N-1,0,1)
# dYdt = function(x){sum((M0$.Ydot(c(0,x),M0$Theta[-M0$i_init]))^2)}
# Y_eq = .fit(Y_eq,dYdt,100)
# Y_eq_scaled = (Y_eq$par-(-2))/(2-(-2))
# points(Y_eq_scaled[1],Y_eq_scaled[2],col="red",pch=8,cex=2)

dev.off()

#
###
