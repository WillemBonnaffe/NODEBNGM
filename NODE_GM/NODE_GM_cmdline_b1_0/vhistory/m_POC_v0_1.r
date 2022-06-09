###########
## m_POC ##
###########

## goal: perform assessment of quality of inference of NODE by gradient matching 

## update:
## 03-03-2022 - created v0_0 from m_NODE_GM_v0_17.r
## 07-03-2022 - created v0_1
##            - re-introduced scaling of variables to help with the fitting 

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
N_e    	 = 3      	       # number of ensemble elements
W_o    	 = 10    	       # number of nodes in observation model
W_p    	 = 10     	       # number of nodes in process model
sd1_o  	 = 0.01   	       # sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    	       # sd in priors of observation model
sd1_p  	 = 0.01   	       # sd in likelihood of residuals of process model
sd2_p  	 = 0.1  	       # sd in priors of process model
alpha_i  = 2      	       # factor of interpolation
alpha_o  = nrow(TS)-1      # scaling factor for observation model

## initiate NODE
t     = TS[,1]                                    # vector of time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)    # normalise interpolated time steps
Y     = as.matrix(TS[,-1])                                   # matrix containing time series
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

## standardise data
Y = log(Y)
mean_ = mean(Y)
sd_ = sd(Y)
Y   = (Y-mean_)/sd_

## fit observation model
Omega_o = list()
for (i in 1:ncol(Y))
{
	Omega_o_      = rnorm(W_o*3,0,sd2_o)
	Omega_o_      = argmax.logPost_o(t,Y[,i],Omega_o_,sd1_o,sd2_o)
	Omega_o[[i]]  = Omega_o_
	message(logPost_o(t,Y[,i],Omega_o_,sd1_o,sd2_o))
}

# ## fit via DEMCpp
# library(Rcpp)
# sourceCpp("cpp/DEMCpp_v0.1.cpp")
# # Omega_o = list()
# chainList = list()
# for(i in 1:ncol(Y))
# {
# 	chainList[[i]] = DEMCpp(list(dTarget=function(x)logPost_o(t,Y[,i],x,sd1_o,sd2_o),
# 						Theta_0=Omega_o[[i]], # Theta_0=rnorm(W_o*3,0,sd2_o),
# 						epsilon=0.0001,
# 						nIt=100000))$chain
# 	Omega_o[[i]] = chainList[[i]][which.max(chain[,1]),-1]
# 	message(logPost_o(t,Y[,i],Omega_o[[i]],sd1_o,sd2_o))
# }
# 
# ## evaluate ensemble
# for(i in 1:ncol(Y))
# {
# 	s        = seq(1,nrow(chainList[[i]]),nrow(chainList[[i]])/100)
# 	tmp      = t(apply(chainList[[i]][s,-1],1,function(x) f_o.eval(t_,x)))
# 	E.Ybar_o = apply(tmp,2,mean)
# 	q05.Ybar_o = apply(tmp,2,quantile,p=0.05)
# 	q95.Ybar_o = apply(tmp,2,quantile,p=0.95)
# 	#
# 	plot(t,Y[,i])
# 	lines(t_,E.Ybar_o)
# 	lines(t_,q05.Ybar_o,lty=2)
# 	lines(t_,q95.Ybar_o,lty=2)
# }

## fit via DEMCO
library(Rcpp)
sourceCpp("cpp/DEMCO_v0.5.cpp")
# Omega_o = list()
for(i in 1:ncol(Y))
{
	chain = DEMCOpp(list(dTarget=function(x)logPost_o(t,Y[,i],x,sd1_o,sd2_o),
						Theta_0=Omega_o[[i]], # Theta_0=rnorm(W_o*3,0,sd2_o),
						epsilon=0.0001,
						lambda=100,
						nIt=10000))$chainList
	Omega_o[[i]] = chain[which.max(chain[,1]),-1]
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

## visualise periodic components
for( i in 1:ncol(Y))
{
	plot(t_,Ybar_o[,i],type="l")
	Ybar_o_ = t(apply(matrix(Omega_o[[i]],ncol=3),1,function(x) f_o.eval(t_,x)))
	for(j in 1:W_o)
	{
		lines(t_,Ybar_o_[j,],col=rainbow(W_o)[j]) 
	}
}

## unscale
Y          = exp(mean_ + sd_*Y)
Ybar_o     = exp(mean_ + sd_*Ybar_o)
ddt.Ybar_o = sd_*Ybar_o*ddt.Ybar_o 

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

# #######################
# ## FIT PROCESS MODEL ##
# #######################
# 
# ## standardise Ybar_o
# mean_ = mean(Ybar_o)
# sd_   = sd(Ybar_o)
# std.Ybar_o = (Ybar_o-mean_)/sd_
# 
# ## standardise Ybar_o
# mean_ = mean(Ybar_true)
# sd_   = sd(Ybar_true)
# std.Ybar_true = (Ybar_true-mean_)/sd_
# 
# 
# ## fit observation model
# Omega_p = list()
# for (i in 1:ncol(Y))
# {
# 	Omega_p_      = rnorm(W_p*(2+ncol(Y)),0,sd2_p)
# 	Omega_p_      = argmax.logPost_p(std.Ybar_true,1/Ybar_true[,i]*ddt.Ybar_true[,i],Omega_p_,sd1_p,sd2_p)
# 	Omega_p[[i]]  = Omega_p_
# }
# 
# ## evaluate observation model
# ddt.Ybar_p     = NULL
# ddx.ddt.Ybar_p = NULL
# for(i in 1:ncol(Y))
# {
# 	ddt.Ybar_p_     = ddt.f_p.eval(std.Ybar_true,Omega_p[[i]])
# 	ddt.Ybar_p      = cbind(ddt.Ybar_p,ddt.Ybar_p_)
# 	ddx.ddt.Ybar_p_ = 1/sd_*t(ddx.ddt.f_p.eval(std.Ybar_true,Omega_p[[i]]))
# 	ddx.ddt.Ybar_p  = cbind(ddx.ddt.Ybar_p,ddx.ddt.Ybar_p_)
# }
# 
# ## visualise
# par(mfrow=c(2,1))
# for(i in 1:ncol(Y))
# {
# 	##
# 	plot(t_true,1/Ybar_true[,i]*ddt.Ybar_true[,i],lty=2,type="l")
# 	lines(t_,1/Ybar_o[,i]*ddt.Ybar_o[,i],type="l",col="black")
# 	lines(t_true,ddt.Ybar_p[,i],type="l",col="red")
# 
# 	##
# 	plot(t_true,ddx.ddt.Ybar_p[,1+(i-1)*2],ylim=c(-2,2),col="white")
# 	for(j in 1:ncol(Y))
# 	{
# 		lines(t_true,ddx.ddt.Ybar_p[,j+(i-1)*2],type="l",col=j)
# 		lines(t_true,ddx.r_true[,j+(i-1)*2],type="l",col=j,lty=2)
# 	}	
# }
# par(mfrow=c(1,1))
# 
# #
# ###

dev.off()

###########################
## FIT OBSERVATION MODEL ##
###########################

## method 1 ##
## source: https://stackoverflow.com/questions/41435777/perform-fourier-analysis-to-a-time-series-in-r

# nff = function(x = NULL, n = NULL, up = 10L, plot = TRUE, add = FALSE, main = NULL, ...){
#   #The direct transformation
#   #The first frequency is DC, the rest are duplicated
#   dff = fft(x)
#   #The time
#   t = seq(from = 1, to = length(x))
#   #Upsampled time
#   nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
#   #New spectrum
#   ndff = array(data = 0, dim = c(length(nt), 1L))
#   ndff[1] = dff[1] #Always, it's the DC component
#   if(n != 0){
#     ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
#     #The negative ones are trickier
#     ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
#   }
#   #The inverses
#   indff = fft(ndff/length(t), inverse = TRUE)
#   idff = fft(dff/length(t), inverse = TRUE)
#   if(plot){
#     if(!add){
#       plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
#         main = ifelse(is.null(main), paste(n, "harmonics"), main))
#       lines(y = Mod(idff), x = t, col = adjustcolor(1L, alpha = 0.5))
#     }
#     lines(y = Mod(indff), x = nt, ...)
#   }
#   ret = data.frame(time = nt, y = Mod(indff))
#   return(ret)
# }

# par(mfrow=c(2,3))
# #
# y     <- Y[,1] 
# rg    <- diff(range(y))
# 
# res = nff(x = y, n = 18L, up = 100L, col = 2L)
# 
# colors = rainbow(36L, alpha = 0.3)
# nff(x = y, n = 36L, up = 100L, col = colors[1])
# for(i in 1:18){
#   ad = ifelse(i == 1, FALSE, TRUE)
#   nff(x = y, n = i, up = 100L, col = colors[i], add = ad, main = "All waves up to 18th harmonic")
# }
# 
# # t     <- 1:73
# y     <- Y[,2] 
# rg    <- diff(range(y))
# 
# res = nff(x = y, n = 18L, up = 100L, col = 2L)
# 
# colors = rainbow(36L, alpha = 0.3)
# nff(x = y, n = 36L, up = 100L, col = colors[1])
# for(i in 1:18){
#   ad = ifelse(i == 1, FALSE, TRUE)
#   nff(x = y, n = i, up = 100L, col = colors[i], add = ad, main = "All waves up to 18th harmonic")
# }
# #
# mfrow=c(1,1)

# ## method 2 ## 
# ## source: http://www.di.fc.ul.pt/~jpn/r/fourier/fourier2.html
# #
# library(pracma) # sech, tanh, fftshift, ifft
# #
# L <- 20;  # interval of interest
# n <- 40;  # sampling rate
# # 
# x = head(seq(0,L,len=n+1), -1);  # the sampling rate was taken on regular intervals
# k <- (2*pi/L) * c(seq(0,n/2-1),seq(-n/2,-1)); # rescaling frequencies from [0,2pi] to [-L/2,L/2]
# # 
# Ybar_true      = ode(y=Y_true_0,times=x,func=ddt.Y_true,parms=p)
# ddt.Ybar_true  = t(apply(Ybar_true[,2:3],1,function(x)unlist(ddt.Y_true(NULL,x,p))))
# #
# u    = nff(Ybar_true[,3],n=5L)[,2] # a function
# t    = nff(Ybar_true[,3],n=5L)[,1]
# ud   = ddt.Ybar_true[,1]  # its first derivative
# #
# ut   = fft(u)
# uds  = ifft( (1i*k)  *ut)     # ifft: Inverse fast Fourier transform
# #
# # draw first derivative and the fourier approximation
# par(mfrow=c(2,1))
# plot(t,u,type='l',col='black', xlab="x", ylab="f(x)")
# lines(t,ifft(ut),col="red")
# plot(x*2.05+1,ud,type='l',col='black', xlab="x", ylab="df/dx")
# lines(t,uds,col="red")
# par(mfrow=c(1,1))

# ## method 2 ## 
# ## source: http://www.di.fc.ul.pt/~jpn/r/fourier/fourier2.html
# #
# library(pracma) # sech, tanh, fftshift, ifft
# #
# L <- 20;  # interval of interest
# n <- 40;  # sampling rate
# # 
# x = head(seq(0,L,len=n+1), -1);  # the sampling rate was taken on regular intervals
# k <- (2*pi/L) * c(seq(0,n/2-1),seq(-n/2,-1)); # rescaling frequencies from [0,2pi] to [-L/2,L/2]
# # 
# Ybar_true      = ode(y=Y_true_0,times=x,func=ddt.Y_true,parms=p)
# ddt.Ybar_true  = t(apply(Ybar_true[,2:3],1,function(x)unlist(ddt.Y_true(NULL,x,p))))
# u    = Ybar_true[,2] # a function
# ud   = ddt.Ybar_true[,1]  # its first derivative
# #
# ut   = fft(u)
# uds  = ifft( (1i*k)  *ut)     # ifft: Inverse fast Fourier transform
# #
# # draw first derivative and the fourier approximation
# par(mfrow=c(2,1))
# plot(x,u,type='l',col='black', xlab="x", ylab="f(x)")
# lines(x,ifft(ut),col="red")
# plot(x,ud,type='l',col='black', xlab="x", ylab="df/dx")
# lines(x,uds,col="red")
# par(mfrow=c(1,1))

#
###

