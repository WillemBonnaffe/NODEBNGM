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
## 17-02-2022 - created v0_13
##            - cleaned code
##            - improved visualisations
## 18-02-2022 - created v0_14
##            - changed fitting function
##            - code summary function (do not calculate quantities on the fly)
##            - separate functions for fitting the observation and process model

## next:
## - code sample.NODEGM
## - code predict.NODEGM
## - code storage of chain/ensemble
## - cool visualisations
## - separate fitting script from visualisation script
## - function init (?)
## - clean code


##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_33.r")

## load data
TS = read.table("TS.csv",sep=";",header=T)

## transform positive variables
TS[,c("Hare","Lynx")] = log(TS[,c("Hare","Lynx")])

## parameters
N_e    	 = 3      # number of ensemble elements
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

###########
## CHAIN ##
###########

## goal: 

## chain
chain = list()
for(k in 1:N_e)
{
	message(paste("Sample ",k,"/",N_e,sep=""))
	chain[[k]] = fit.NODEGM(t,Y,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p,alpha_i,alpha_o) 
}

##
pred = lapply(chain,function(x) eval.NODEGM(t_,x))

##
results = summary.NODEGM(pred,t_,N)
attach(results)

#
###

####################
## VISUALISATIONS ##
####################

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

##
Ybar_o = E.Ybar_o
ddt.Ybar_o = E.ddt.Ybar_o
ddt.Ybar_p = E.ddt.Ybar_p
ddx.ddt.Ybar_p = E.ddx.ddt.Ybar_p
Geber_p        = E.Geber_p

pdf("results_2.pdf")
#
## visualise model
for(i in 1:N)
{
	par(mfrow=c(2,2))
	#
	## interpolation
	plot(t_,rep(0,length(t_)),ylim=c(-2,2),type="l",lty=2)
	points(t,Y[,i])
	lines(t_,Ybar_o[,i],col=rainbow(N)[i])
	#
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	lines(t_,ddt.Ybar_o[,i],col="black")
	lines(t_,ddt.Ybar_p[,i],col=rainbow(N)[i])
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	for(j in 1:N) lines(t_, ddx.ddt.Ybar_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1),type="l",lty=2)
	for(j in 1:N) lines(t_, Geber_p[,(i-1)*N+(1:N)][,j],col=rainbow(N)[j])
	#
	par(mfrow=c(1,1))
}
#
##
par(mfrow=c(2,2))
#
for(i in 1:N)
{
	for(j in 1:N)
	{
		for(k in 1:N) 
			{
				plot(Ybar_o[,k],ddx.ddt.Ybar_p[,(i-1)*N + 1:N][,j],ylim=c(-1,1),main=paste(j,"->",i," : f(",k,")",sep=""))
				lines(c(-1,1)*10,c(0,0),lty=2)
			}
	}
}
#
par(mfrow=c(1,1))
#
par(mfrow=c(N,N))
#
for(i in 1:N)
{
	for(j in 1:N) 
		{
			x = Geber_p[,(i-1)*N + 1:N]
			y = Geber_p[,(j-1)*N + 1:N]
			plot(x,y,type="l")
			lines(c(-1,1)*10,c(0,0),lty=2)
			lines(c(0,0),c(-1,1)*10,lty=2)
		}
}
#
par(mfrow=c(1,1))
#
dev.off()

#
###
