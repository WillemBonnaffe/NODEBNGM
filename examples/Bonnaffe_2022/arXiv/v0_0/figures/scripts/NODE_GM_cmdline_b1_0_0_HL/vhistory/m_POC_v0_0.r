###########
## m_POC ##
###########

## goal: perform assessment of quality of inference of NODE by gradient matching 

## update:
# 03-03-2022 - created v0_0 from m_NODE_GM_v0_17.r

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
p          = c(1.0,0.6,0.1,1.0)
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

## visualise
plot(TS[,1],TS[,2],col="red",type="l")
lines(TS[,1],TS[,3],col="blue")

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
N_e    	 = 3      	       # number of ensemble elements
W_o    	 = 50    	       # number of nodes in observation model
W_p    	 = 10     	       # number of nodes in process model
sd1_o  	 = 0.01   	       # sd in likelihood of residuals of observation model
sd2_o  	 = 0.01    	       # sd in priors of observation model
sd1_p  	 = 0.01   	       # sd in likelihood of residuals of process model
sd2_p  	 = 0.1  	       # sd in priors of process model
alpha_i  = 2      	       # factor of interpolation
alpha_o  = nrow(TS)-1      # scaling factor for observation model

## initiate NODE
t     = TS[,1]                                    # vector of time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)    # normalise interpolated time steps
Y     = TS[,-1]                                   # matrix containing time series
#
# ## normalise time vector
# std.t  = (t-min(t))/(max(t)-min(t))*alpha_o                         # normalise time steps
# std.t_ = seq(min(std.t),max(std.t),(std.t[2]-std.t[1])/alpha_i)     # normalise interpolated time steps
#
# ## standardise state variables
# mean_ = apply(Y,2,mean)                         # save mean of each variable in time series
# sd_   = apply(Y,2,sd)                           # save sd of each variable in time series
# std.Y = t(apply(Y,1,function(x)(x-mean_)/sd_))  # standardise each variable
#
# mean_ = mean(as.matrix(Y))
# sd_   = sd(as.matrix(Y))
# std.Y = (Y-mean_)/(sd_)

#
###

pdf("results.pdf",height=16,width=16)

###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model

## fit observation model
Omega_o = list()
for (i in 1:ncol(Y))
{
	Omega_o_      = rnorm(W_o*3,0,sd2_o)
	Omega_o_      = argmax.logPost_o(t,Y[,i],Omega_o_,sd1_o,sd2_o)
	Omega_o[[i]]  = Omega_o_
}

## evaluate observation model
Ybar_o     = NULL 
ddt.Ybar_o = NULL
for(i in 1:ncol(Y))
{
	Ybar_o_     = f_o.eval(t_,Omega_o[[i]])
	Ybar_o      = cbind(Ybar_o,Ybar_o_)
	ddt.Ybar_o_ = ddt.f_o.eval(t_,Omega_o[[i]])
	ddt.Ybar_o  = cbind(ddt.Ybar_o,ddt.Ybar_o_)
}

## visualise
par(mfrow=c(2,1))
for(i in 1:ncol(Y))
{
	##
	plot(t,Y[,i])
	lines(t_,Ybar_o[,i],col="red")
	lines(t_true,Ybar_true[,i],lty=2)
	#
	##
	plot(t_true,ddt.Ybar_true[,i],lty=2,type="l")
	lines(t_,ddt.Ybar_o[,i],type="l",col="red")
}
par(mfrow=c(1,1))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## standardise Ybar_o
mean_ = mean(Ybar_o)
sd_   = sd(Ybar_o)
std.Ybar_o = (Ybar_o-mean_)/sd_

## standardise Ybar_o
mean_ = mean(Ybar_true)
sd_   = sd(Ybar_true)
std.Ybar_true = (Ybar_true-mean_)/sd_


## fit observation model
Omega_p = list()
for (i in 1:ncol(Y))
{
	Omega_p_      = rnorm(W_p*(2+ncol(Y)),0,sd2_p)
	Omega_p_      = argmax.logPost_p(std.Ybar_true,1/Ybar_true[,i]*ddt.Ybar_true[,i],Omega_p_,sd1_p,sd2_p)
	Omega_p[[i]]  = Omega_p_
}

## evaluate observation model
ddt.Ybar_p     = NULL
ddx.ddt.Ybar_p = NULL
for(i in 1:ncol(Y))
{
	ddt.Ybar_p_     = ddt.f_p.eval(std.Ybar_true,Omega_p[[i]])
	ddt.Ybar_p      = cbind(ddt.Ybar_p,ddt.Ybar_p_)
	ddx.ddt.Ybar_p_ = 1/sd_*t(ddx.ddt.f_p.eval(std.Ybar_true,Omega_p[[i]]))
	ddx.ddt.Ybar_p  = cbind(ddx.ddt.Ybar_p,ddx.ddt.Ybar_p_)
}

## visualise
par(mfrow=c(2,1))
for(i in 1:ncol(Y))
{
	##
	plot(t_true,1/Ybar_true[,i]*ddt.Ybar_true[,i],lty=2,type="l")
	lines(t_,1/Ybar_o[,i]*ddt.Ybar_o[,i],type="l",col="black")
	lines(t_true,ddt.Ybar_p[,i],type="l",col="red")

	##
	plot(t_true,ddx.ddt.Ybar_p[,1+(i-1)*2],ylim=c(-2,2),col="white")
	for(j in 1:ncol(Y))
	{
		lines(t_true,ddx.ddt.Ybar_p[,j+(i-1)*2],type="l",col=j)
		lines(t_true,ddx.r_true[,j+(i-1)*2],type="l",col=j,lty=2)
	}	
}
par(mfrow=c(1,1))

#
###

dev.off()
