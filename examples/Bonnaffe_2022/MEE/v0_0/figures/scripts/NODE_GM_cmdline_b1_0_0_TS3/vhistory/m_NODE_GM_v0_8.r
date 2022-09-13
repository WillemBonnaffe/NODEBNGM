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

## next:
## test the scaling of the covariates to make the marginal likelihood appraoch work for the process model

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
N_e    	 = 3      # number of ensemble elements
W_o    	 = 100    # number of nodes in observation model
W_p    	 = 100    # number of nodes in process model
sd1_o  	 = 0.01   # sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    # sd in priors of observation model
sd1_p  	 = 0.01   # sd in likelihood of residuals of process model
sd2_p  	 = 0.1    # sd in priors of process model
alpha_o  = 100    # scaling factor for observation model
alpha_p  = 100    # scaling factor for process model
method   = 1      # 1 for normal posterior / 2 for marginal posterior

#
###

##########
## MAIN ##
##########

## goal: fit the NODE model

## plot time series
plot(TS[,1],TS[,2],xlab="Time (days)",ylab="Density (normalised)")
for(j in 1:2)
{
	points(TS[,1],TS[,j+1])
	lines(TS[,1],TS[,j+1],col=rainbow(ncol(TS)-1)[j])
}

## fit NODE
fit.NODE(TS,N_e,W_o,W_p,sd1_o,sd2_o,sd1_p,sd2_p,alpha_o,alpha_p,method)

#
###
