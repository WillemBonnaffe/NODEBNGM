############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

## update:
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
#            - removed unused sections
# 22-11-2021 - crated v0_7
#            - transferred main functions to function file
# 23-11-2021 - created v0_8
#            - implemented scaling of the process model
# 26-11-2021 - created v0_9
#            - implemented a fitting method that does not require anchor sampling
# 29-11-2021 - created v0_10
#            - experimented further with code

## next:
## re-implement the Geber method
## clean code
## cool visualisations

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_30.r")

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
sd2_p  	 = 1.0    # sd in priors of process model
alpha_i  = 1      # factor of interpolation
alpha_o  = 100    # scaling factor for observation model
alpha_p  = 1      # scaling factor for process model

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
# # ## plot time series
# # plot(TS[,1],TS[,2],xlab="Time (days)",ylab="Density (normalised)")
# # for(j in 1:2)
# # {
# # 	points(TS[,1],TS[,j+1])
# # 	lines(TS[,1],TS[,j+1],col=rainbow(ncol(TS)-1)[j])
# # }
# 
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
# Ybar_o = Ybar_o * alpha_p
# ddx.ddt.Ybar_p = ddt.Ybar_p = NULL
# for (i in 1:N)
# {
# 	Omega_p = rnorm(W_p*(2 + N),0,sd2_p)
# 	Omega_p = argmax.logPost_p(Ybar_o,ddt.Ybar_o[,i],Omega_p)
# 	ddt.Ybar_p = cbind(ddt.Ybar_p,ddt.f_p.eval(Ybar_o,Omega_p))
# 	ddx.ddt.Ybar_p = cbind(ddx.ddt.Ybar_p,t(ddx.ddt.f_p.eval(Ybar_o,Omega_p)))
# }
# 
# pdf("results.pdf")
# #
# ## visualise model
# for(i in 1:N)
# {
# 	par(mfrow=c(3,1))
# 	#
# 	## interpolation
# 	plot(t_,rep(0,length(t_)),ylim=c(-2,2),type="l",lty=2)
# 	points(t,Y[,i])
# 	lines(t_,Ybar_o[,i]/alpha_p,col=rainbow(N)[i])
# 	#
# 	## dynamics
# 	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
# 	lines(t_,ddt.Ybar_o[,i],col="black")
# 	lines(t_,ddt.Ybar_p[,i],col=rainbow(N)[i])
# 	#
# 	## effects
# 	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
# 	for(j in 1:N) lines(t_, ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j] * alpha_p,col=rainbow(N)[j])
# 	#
# 	par(mfrow=c(1,1))
# }
# #
# ##
# for(i in 1:N) for(j in 1:N) for(k in 1:N) plot(Ybar_o[,i],ddx.ddt.Ybar_p[,(j-1)*N + 1:N][,k]*alpha_p)
# #
# dev.off()
# 
# 
# ## fit NODE
# # fit.NODE(TS,N_e,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p,alpha_o,alpha_p,method)
# 
# #
# ###

###########
## CHAIN ##
###########

## goal: 

## initiate chain
chain = list()

for(k in 1:N_e)
{
	message(paste("Sample ",k,"/",N_e,sep=""))

	## observation model
	Ybar_o = ddt.Ybar_o = NULL
	for (i in 1:N)
	{
		Omega_o = rnorm(W_o*3,0,sd2_o)
		Omega_o = argmax.logPost_o(t,Y[,i],Omega_o)
		Ybar_o = cbind(Ybar_o,f_o.eval(t_,Omega_o))	
		ddt.Ybar_o = cbind(ddt.Ybar_o,ddt.f_o.eval(t_,Omega_o))	
	}
	
	## process model
	Ybar_o = Ybar_o * alpha_p
	ddx.ddt.Ybar_p = ddt.Ybar_p = NULL
	for (i in 1:N)
	{
		Omega_p = rnorm(W_p*(2 + N),0,sd2_p)
		Omega_p = argmax.logPost_p(Ybar_o,ddt.Ybar_o[,i],Omega_p)
		ddt.Ybar_p = cbind(ddt.Ybar_p,ddt.f_p.eval(Ybar_o,Omega_p))
		ddx.ddt.Ybar_p = cbind(ddx.ddt.Ybar_p,t(ddx.ddt.f_p.eval(Ybar_o,Omega_p)))
	}

	## store chain element
	chain[[k]] = cbind(Ybar_o,ddt.Ybar_o,ddt.Ybar_p,ddx.ddt.Ybar_p)
}

## compute mean and error around predictions
chain_ = matrix(unlist(chain),ncol=N_e)
E.chain   = matrix(apply(chain_,1,mean),nrow=n_)
q05.chain = matrix(apply(chain_,1,quantile,p=0.05),nrow=n_) 
q95.chain = matrix(apply(chain_,1,quantile,p=0.95),nrow=n_)

## compute expectation of variables
E.Ybar_o = E.chain[,1:N]
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

pdf("results.pdf")
#
## visualise model
for(i in 1:N)
{
	par(mfrow=c(3,1))
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
	for(j in 1:N) lines(t_, E.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j] * alpha_p,col=rainbow(N)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],rev(q95.ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	#
	par(mfrow=c(1,1))
}
#
dev.off()


## fit NODE
# fit.NODE(TS,N_e,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p,alpha_o,alpha_p,method)

#
###
