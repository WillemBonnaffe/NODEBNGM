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
## 24-02-2022 - created v0_16
##            - cleaned code
##            - code storage of chain/ensemble
##            - clean code (esp. dependences, limit dependence on global variables)
## 03-03-2022 - created v0_17
##            - simplified preparation of data
## 16-03-2022 - created v0_18
##            - separate fitting script from visualisation script
##            - update fitting and plotting function with new repository
## 17-03-2022 - created v0_19
##            - impemented relative contributions
##            - implemented diamond plot
## 18-03-2022 - extra visualisations (diamond plots, phase space, ...)
##            - improve figures
## 20-03-2022 - created v0_20
##            - cleaned code by separating observation and process model

## next:
##            - allow for different values of priors for each time series

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_47.r")

## load data
# TS = as.matrix(read.table("data/TS_LV_2.csv",sep=";",header=T))
#
TS = read.table("data/TS_RPS.csv",sep=";",header=T)
s  = seq(0,600,10)
TS = TS[s,]
TS[,1] = TS[,1]/10
TS[,-1] = TS[,-1]/1000

# TS = read.table("data/TS_HL.csv",sep=";",header=T)
# TS[,1] = (TS[,1]-min(TS[,1]))
# TS[,-1] = TS[,-1]/10

## initiate time series 
t     = TS[,1]                                    # vector of time steps
Y     = TS[,-1]                                   # matrix containing time series

## output file
pdf("results.pdf")

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model 

## parameters
N_o    	 = 3      	       # number of ensemble elements
W_o    	 = 10    	       # number of nodes in observation model
sd1_o  	 = 0.1   	       # sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    	       # sd in priors of observation model
alpha_i  = 1     	       # factor of interpolation

## interpolated time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)

## fit model_o 
model_o = fit.NODEGM_o(t,t_,Y,W_o,sd1_o,sd2_o,N_o,log=T,logMar=T)
# model_o = DFT.NODEGM_o(t,t_,Y,W_o,log=T)

## visualise
plot.NODEGM_o(t,t_,Y,model_o,col=rainbow(ncol(Y)))

## prepare data
Yhat_o     = NULL
ddt.Yhat_o = NULL
for(i in 1:(ncol(TS)-1))
{
    Yhat_o     = cbind(Yhat_o,    model_o[[i]]$MaP.Yhat_o)
    ddt.Yhat_o = cbind(ddt.Yhat_o,model_o[[i]]$MaP.ddt.Yhat_o)
}

## compute per-capita growth rate
r_o = 1/Yhat_o*ddt.Yhat_o

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model
 
## parameters
N_p    	 = 10      	       # number of ensemble elements
W_p    	 = 10     	       # number of nodes in process model
sd1_p  	 = 0.1   	       # sd in likelihood of residuals of process model
sd2_p  	 = 0.25  	       # sd in priors of process model

## fit model_p
model_p = fit.NODEGM_p(Yhat_o,ddt.Yhat_o,r_o,W_p,sd1_p,sd2_p,N_p)

## visualise results
plot.NODEGM_p(t_,r_o,model_p,col=rainbow(ncol(Y)))
# diamondPlot.NODEGM_p(model_p,labels=colnames(Y))

## terminate
# dev.off()

#
###

# ########################
# ## DEMC PROCESS MODEL ##
# ########################
# 
# ## goal: fit process model with DEMC
# 
# ## parameters
# N_p    	 = 10000      	   # number of ensemble elements
# W_p    	 = 1     	       # number of nodes in process model
# sd1_p  	 = 0.1   	       # sd in likelihood of residuals of process model
# sd2_p  	 = 0.3  	       # sd in priors of process model
# 
# ## fit model_p
# par(mfrow=c(3,1))
# #
# for(i in 1:ncol(Y))
# {
# 	model_p   = sample.model_p(Yhat_o,r_o[,i],W_p,sd1_p,sd2_p,N_p)
# 	s = seq(N_p/2,N_p,(N_p-N_p/2)/1000)
# 	model_p$ensemble = model_p$ensemble[s,]
# 	results_p = summary.model_p(Yhat_o,ddt.Yhat_o,model_p)
# 	
# 	## visualise
# 	plot.model_p(t_,r_o[,i],results_p)
# }
# #
# par(mfrow=c(1,1))
# 
# ## terminate
# dev.off()
# 
# #
# ###
