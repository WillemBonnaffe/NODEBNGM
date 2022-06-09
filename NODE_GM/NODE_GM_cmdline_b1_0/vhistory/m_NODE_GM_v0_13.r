############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

## update:
## 10-03-2021 - created version 0.0
## 11-03-2021 - created version 0.1
##            - implemented figures for the paper 
## 29-03-2021 - created version 0.2
##            - combined figures for the paper
## 25-05-2021 - created version 0.3
##            - modularised code to work with any single time series file
## 18-11-2021 - crated version 0_5
##            - cleaned code
## 22-11-2021 - cleaned code
## 22-11-2021 - created v0_6
##            - removed unused sections
## 22-11-2021 - crated v0_7
##            - transferred main functions to function file
## 23-11-2021 - created v0_8
##            - implemented scaling of the process model
## 26-11-2021 - created v0_9
##            - implemented a fitting method that does not require anchor sampling
## 29-11-2021 - created v0_10
##            - experimented further with code
## 30-11-2021 - created v0_11
##            - implemented Geber method
##            - created v0_12
##            - implemented main function
## 18-02-2022 - created v0_13
##            - cleaned code
##            - improved visualisations

## next:
## - code summary function (do not calculate quantities on the fly)
## - clean code
## - cool visualisations
## - code ensemble.NODEGM
## - think about how much detail you want in code (or putting stuff in functions)
## - think about whether to pre-compute derived quantities or to get the ensemble and then compute quantities (the latter, like a summary function)
## - separate fitting the observation and process model again (?)
## - code storage of chain/ensemble
## - separate fitting script from visualisation script


##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_32.r")

## load data
TS = read.table("TS.csv",sep=";",header=T)

## transform positive variables
TS[,c("Hare","Lynx")] = log(TS[,c("Hare","Lynx")])

## parameters
N_e    	 = 10     # number of ensemble elements
W_o    	 = 100    # number of nodes in observation model
W_p    	 = 100    # number of nodes in process model
sd1_o  	 = 0.01   # sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    # sd in priors of observation model
sd1_p  	 = 0.01   # sd in likelihood of residuals of process model
sd2_p  	 = 0.5    # sd in priors of process model
alpha_i  = 1      # factor of interpolation
alpha_o  = 100    # scaling factor for observation model

## initiate NODE
n     = nrow(TS)                                  # number of time steps
n_    = n*alpha_i                                 # number of interpolated time steps
N     = ncol(TS) - 1                              # number of variables
t     = TS[,1]                                    # vector of time steps
Y     = TS[,-1]                                   # matrix containing time series
t     = (t-min(t))/(max(t)-min(t))*alpha_o        # normalise time steps
t_    = seq(min(t),max(t),(max(t)-min(t))/n_)     # normalise interpolated time steps
n_    = length(t_)                                # double check number of interpolated time steps
mean_ = apply(Y,2,mean)                           # save mean of each variable in time series
sd_   = apply(Y,2,sd)                             # save sd of each variable in time series
Y     = t(apply(Y,1,function(x)(x-mean_)/sd_))    # standardise each variable

#
###

# ##########
# ## MAIN ##
# ##########
# 
# ## goal: fit the NODE model
#  
# ## observation model
# Ybar_o = ddt.Ybar_o = NULL
# for (i in 1:N)
# {
# 	Omega_o = rnorm(W_o*3,0,sd2_o)
# 	Omega_o = argmax.logPost_o(t,Y[,i],Omega_o)
# 	Ybar_o = cbind(Ybar_o,f_o.eval(t_,Omega_o))	
# 	ddt.Ybar_o = cbind(ddt.Ybar_o,ddt.f_o.eval(t_,Omega_o))	
# }
# 
# ## process model
# ddx.ddt.Ybar_p = ddt.Ybar_p = NULL
# for (i in 1:N)
# {
# 	Omega_p = rnorm(W_p*(2 + N),0,sd2_p)
# 	Omega_p = argmax.logPost_p(Ybar_o,ddt.Ybar_o[,i],Omega_p)
# 	ddt.Ybar_p = cbind(ddt.Ybar_p,ddt.f_p.eval(Ybar_o,Omega_p))
# 	ddx.ddt.Ybar_p = cbind(ddx.ddt.Ybar_p,t(ddx.ddt.f_p.eval(Ybar_o,Omega_p)))
# }
# 
# ## Geber method
# Geber_p = t(apply(cbind(ddt.Ybar_p,ddx.ddt.Ybar_p),1,function(x) x[1:N] * x[-c(1:N)]))
# 
# 
# pdf("results.pdf")
# #
# ## visualise model
# for(i in 1:N)
# {
# 	par(mfrow=c(2,2))
# 	#
# 	## interpolation
# 	plot(t_,rep(0,length(t_)),ylim=c(-2,2),type="l",lty=2)
# 	points(t,Y[,i])
# 	lines(t_,Ybar_o[,i],col=rainbow(N)[i])
# 	#
# 	## dynamics
# 	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
# 	lines(t_,ddt.Ybar_o[,i],col="black")
# 	lines(t_,ddt.Ybar_p[,i],col=rainbow(N)[i])
# 	#
# 	## effects
# 	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
# 	for(j in 1:N) lines(t_, ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
# 	#
# 	## Geber 
# 	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
# 	for(j in 1:N) lines(t_, Geber_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
# 	#
# 	par(mfrow=c(1,1))
# }
# #
# ##
# par(mfrow=c(2,2))
# #
# for(i in 1:N)
# {
# 	for(j in 1:N)
# 	{
# 		for(k in 1:N) 
# 			{
# 				plot(Ybar_o[,k],ddx.ddt.Ybar_p[,(i-1)*N + 1:N][,j],ylim=c(-1,1),main=paste(j,"->",i," : f(",k,")",sep=""))
# 				lines(c(-1,1)*10,c(0,0),lty=2)
# 			}
# 	}
# }
# #
# par(mfrow=c(1,1))
# #
# par(mfrow=c(1,1))
# #
# par(mfrow=c(N,N))
# #
# for(i in 1:N)
# {
# 	for(j in 1:N) 
# 		{
# 			x = Geber_p[,(i-1)*N + 1:N]
# 			y = Geber_p[,(j-1)*N + 1:N]
# 			plot(x,y,type="l")
# 			lines(c(-1,1)*10,c(0,0),lty=2)
# 			lines(c(0,0),c(-1,1)*10,lty=2)
# 		}
# }
# #
# par(mfrow=c(1,1))
# #
# dev.off()
# 
# #
# # ###

###########
## CHAIN ##
###########

## goal: 

# ## chain
# chain = list()
# for(k in 1:N_e)
# {
# 	message(paste("Sample ",k,"/",N_e,sep=""))
# 	chain[[k]] = fit.NODEGM(TS,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p,alpha_i,alpha_o) 
# }

## compute mean and error around predictions
chain_    = matrix(unlist(chain),ncol=N_e)
E.chain   = matrix(apply(chain_,1,mean),nrow=n_)
q05.chain = matrix(apply(chain_,1,quantile,p=0.05),nrow=n_) 
q95.chain = matrix(apply(chain_,1,quantile,p=0.95),nrow=n_)

## compute expectation of variables
E.Ybar_o   = E.chain[,1:N]
q05.Ybar_o = q05.chain[,1:N]
q95.Ybar_o = q95.chain[,1:N]
#
E.ddt.Ybar_o = E.chain[,1:N + N]
q05.ddt.Ybar_o = q05.chain[,1:N+N]
q95.ddt.Ybar_o = q95.chain[,1:N+N]
# 
E.ddt.Ybar_p = E.chain[,1:N + N*2]
q05.ddt.Ybar_p = q05.chain[,1:N+N*2]
q95.ddt.Ybar_p = q95.chain[,1:N+N*2]
#
E.ddx.ddt.Ybar_p = E.chain[,1:(N*N) + N*3]
q05.ddx.ddt.Ybar_p = q05.chain[,1:(N*N)+N*3]
q95.ddx.ddt.Ybar_p = q95.chain[,1:(N*N)+N*3]
#
E.Geber_p = E.chain[,1:(N*N) + N*3 + N^2]
q05.Geber_p = q05.chain[,1:(N*N) + N*3 + N^2]
q95.Geber_p = q95.chain[,1:(N*N) + N*3 + N^2]

pdf("results.pdf")
#
## visualise model
for(i in 1:N)
{
	par(mfrow=c(4,1),mar=c(3,3,1,1))
	#
	## interpolation
	plot(t_,rep(0,length(t_)),ylim=c(-2,2),type="l",lty=2)
	points(t,Y[,i])
	lines(t_,E.Ybar_o[,i],col=rainbow(N)[i])
	polygon(c(t_,rev(t_)),c(q05.Ybar_o[,i],rev(q95.Ybar_o[,i])),col=rainbow(N,alpha=0.2)[i],border=NA)
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	lines(t_,E.ddt.Ybar_o[,i],col="black")
	lines(t_,E.ddt.Ybar_p[,i],col=rainbow(N)[i])
	polygon(c(t_,rev(t_)),c(q05.ddt.Ybar_p[,i],rev(q95.ddt.Ybar_p[,i])),col=rainbow(N,alpha=0.2)[i],border=NA)
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	for(j in 1:N) lines(t_, E.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],rev(q95.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	for(j in 1:N) lines(t_, E.Geber_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,(i-1)*N+(1:N)][,j],rev(q95.Geber_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	#
	par(mfrow=c(1,1))
}
#
dev.off()

#
###
