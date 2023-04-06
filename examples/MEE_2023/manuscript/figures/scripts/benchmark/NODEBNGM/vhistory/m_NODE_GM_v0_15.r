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
## 23-02-2022 - created v0_15
##            - moved sampling and plotting code to function folder (f_NODE_GM_v0_34.r)
##            - code sample.NODEGM
##            - code predict.NODEGM

## next:
##            - code storage of chain/ensemble
##            - separate fitting script from visualisation script (?)
##            - clean code (esp. dependences, limit dependence on global variables)
##            - extra visualisations (diamond plots, phase space, ...)

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_34.r")

## load data
TS = read.table("TS.csv",sep=";",header=T)

## transform positive variables
TS[,c("Hare","Lynx")] = log(TS[,c("Hare","Lynx")])

## parameters
N_e    	 = 10      	# number of ensemble elements
W_o    	 = 100    	# number of nodes in observation model
W_p    	 = 100     	# number of nodes in process model
sd1_o  	 = 0.01   	# sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    	# sd in priors of observation model
sd1_p  	 = 0.01   	# sd in likelihood of residuals of process model
sd2_p  	 = .1    	# sd in priors of process model
alpha_i  = 2      	# factor of interpolation
alpha_o  = 100    	# scaling factor for observation model

## initiate NODE
t     = TS[,1]                                    # vector of time steps
Y     = TS[,-1]                                   # matrix containing time series
#
N     = ncol(Y)                                   # number of variables
n     = length(t)                                 # number of time steps
n_    = n*alpha_i                                 # number of interpolated time steps
#
## normalise time vector
t_real = t 
t      = (t-min(t))/(max(t)-min(t))*alpha_o        # normalise time steps
t_     = seq(min(t),max(t),(max(t)-min(t))/n_)     # normalise interpolated time steps
n_     = length(t_)                                # double check number of interpolated time steps
#
## standardise state variables
# mean_ = apply(Y,2,mean)                         # save mean of each variable in time series
# sd_   = apply(Y,2,sd)                           # save sd of each variable in time series
# Y     = t(apply(Y,1,function(x)(x-mean_)/sd_))  # standardise each variable
#
mean_ = mean(as.matrix(Y))
sd_   = sd(as.matrix(Y))
Y     = (Y-mean_)/(sd_)

#
###

#########
## FIT ##
#########

## goal: 

## 
# chain    = sample.NODEGM(N_e,t,t_,Y,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p) 
ensemble = getEnsemble.NODEGM(t_,chain) 
results  = summary.NODEGM(ensemble,t_,N)
save(results,file="results.RData")

##
pdf("results.pdf")
plot.NODEGM(TS,results)
dev.off()

#
###
