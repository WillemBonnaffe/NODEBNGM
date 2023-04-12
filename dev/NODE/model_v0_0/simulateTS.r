##################
## simulateTS.r ##
##################

## goal: simulate a time series to test the NODE implementation

##############
## INITIATE ##
##############

## load libraries
library(deSolve)

#
###

###################
## SIMULATE DATA ##
###################

## define model 
dYdt = function(times,state,param)
{
	return(list(c((1.0 - 0.1 * state[2]) * state[1], (0.1 * state[1] - 1.0) * state[2])))
}

## run model
t = seq(0,20,1)
Y_0 = c(0.2,0.5)
TS = ode(y=Y_0,times=t,parms=NULL,func=dYdt)
colnames(TS) = c("t","H","L")

## plot model
# plot(data.frame(TS))

## store
write.table(TS,"TS.csv",sep=";",row.names=FALSE,quote=FALSE)

#
###
