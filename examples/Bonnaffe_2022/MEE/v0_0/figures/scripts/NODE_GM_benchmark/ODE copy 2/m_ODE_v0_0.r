############
## main.R ##
############

## libraries
library(deSolve)
library(Rcpp)

## source
source("R/f_HBM_v1_0.r")
source("R/f_DEMC_v1_0.R")
sourceCpp("Rcpp/f_DEMCpp_v1_0.cpp")
sourceCpp("Rcpp/f_DEMCO_v1_0.cpp")

## make out directory
pathToOut = "out/"
system(paste("mkdir",pathToOut))

#
###

######################
## LOAD TIME SERIES ##
######################

## goal: load time series

## load data
TS = read.table("data/TS_3DLV.csv",sep=";",header=T)
TS = TS[20:50,]
TS[,1] = TS[,1]-min(TS[,1])
Y  = as.matrix(TS)

#
###

##################
## INITIATE ODE ##
##################

## LV system
dY = function(t, stateVect, paramVect)
{
  return(list(c(
    stateVect[1]*(paramVect[1] + paramVect[2]*stateVect[2] + paramVect[3]*stateVect[3]),
    stateVect[2]*(paramVect[4] + paramVect[5]*stateVect[1] + paramVect[6]*stateVect[3]),
    stateVect[3]*(paramVect[7] + paramVect[8]*stateVect[1] + paramVect[9]*stateVect[2])
  )))
}

## posterior
dLogPost = function(Y, Theta)
{
  ## predictive function
  Ybar = ode(y = exp(Theta[1:3]), times = Y[,1], func=dY, parms=Theta[-c(1:3)], method="ode23")
  
  ## calculate posterior
  return(
    sum(log(dgamma(Y[,2],shape=Ybar[,2],rate=1)) + 
        log(dgamma(Y[,3],shape=Ybar[,3],rate=1)) +
        log(dgamma(Y[,4],shape=Ybar[,4],rate=1))
        )
    + sum(log(dunif(Theta,-10,10)))
  )
}

## initiate
Theta_0 = c(log(Y[1,-1]),runif(9,-.1,.1))
# names(Theta_0) = c("V1","V2","V3","V4")
dLogPost(Y,Theta_0)

## initiate chains
dTarget = function(paramVect)
{
  y = dLogPost(Y=Y, paramVect)
  if(!is.nan(y) & !is.na(y)){return(y)}else{return(-Inf)}
}
dTarget(Theta_0) # checks out

## initiate function
initiate = function()
{
  check = F
  while(check == F)
  {
    Theta_0 = c(log(Y[1,-1]),runif(9,-1,1))
    # names(Theta_0) = c("V1","V2","V3","V4")
    if(dTarget(Theta_0) > -Inf){check = T}
  }
  return(Theta_0)
}
dTarget(initiate())

#
###

################
## ODE CHAINS ##
################

## avail if need to re-run the chains
chainList = list()
timeVect = c(0)

# ## R implementation of EDEMC
# timeVect[1] = system.time(
#     for(i in 1:1)
#     {
#         chainList[[i]] = DEMC(dTarget = dTarget,
#                               Theta_0 = initiate(), # initiate from the priors
#                               gamma   = 2.38/sqrt(2*12),
#                               epsilon = 0.001,
#                               nIt     = 1000,
#                               msg     = T)[-1,]
#     })[3]

## RCpp implementation of EDEMC
timeVect[1] = system.time(
  for(i in 1:1)
  {
    Theta_0 = initiate()
    # Theta_0 = MaP
    chainList[[i]] = DEMCOpp(list(dTarget = dTarget,
                                 Theta_0 = Theta_0, # initiate from the priors
                                 gamma   = 2.38/sqrt(2*12),
                                 epsilon = 0.001,
                                 lambda = 100,
                                 nIt     = 2000))[["chainList"]]
    MaP = chainList.summaryTab(chainList.thin(chainList.burn(chainList,1:500)))$estimates[-1,"MaP"]
  })[3]


## RCpp implementation of EDEMC
timeVect[1] = system.time(
    for(i in 1:1)
    {
        # Theta_0 = initiate()
        Theta_0 = MaP
        chainList[[i]] = DEMCpp(list(dTarget = dTarget,
                                     Theta_0 = Theta_0, # initiate from the priors
                                     gamma   = 2.38/sqrt(2*12),
                                     epsilon = 0.001,
                                     nIt     = 1000))[["chainList"]]
        MaP = chainList.summaryTab(chainList.thin(chainList.burn(chainList,1:500)))$estimates[-1,"MaP"]
    })[3]

## store chains
chainList.write(chainList,"out_ODE")

#
###

################
## COMPARISON ##
################

## load chains
chainList = chainList.read("out_ODE")
# timeVect = c(312.362,125.786,102.705)

## untransform
# for(l in 1:length(chainList)){chainList[[l]][,-1] = exp(chainList[[l]][,-1])}

## traces MC
chainList.tracePlot(chainList)
chainList.tracePlot(chainList.thin(chainList.burn(chainList,1:500)))

## plot posterior distribution DEMC
chainList.postPlot(chainList.thin(chainList.burn(chainList,1:500)),1000)

## ac plots MC
par(mfrow=c(3,4))
chainList.acPlot(chainList.burn(chainList,1:500))
# chainList.acPlot(chainList.thin(chainList.burn(chainList[1:3],1:1000)))
par(mfrow=c(1,1))

## summary table
chainList.summaryTab(chainList.thin(chainList.burn(chainList,1:500)),returnSignif=F)[["estimates"]][-1,]
#
ESSTab = matrix(unlist(chainList.ESS(chainList.thin(chainList.burn(chainList,1:500)))),byrow=T,ncol=4)
ESSVect = round(c(apply(ESSTab,2,mean)),2)
#
# write.table(comparisonTab,"out_ODE/comparisonTab.txt",sep=" & ",quote=F, row.names=F)

## model fit
#
## predictions
MaP = chainList.summaryTab(chainList.thin(chainList.burn(chainList,1:500)))$estimates[-1,"MaP"]
Ybar = ode(y = exp(MaP[1:3]), times = seq(0,max(Y[,1]),0.01), func=dY, parms=MaP[-c(1:3)], method="ode45")
#
## plot
par(mfrow=c(1,1),mar=c(5,5,0,0),oma=c(0,0,1,1))
k = 1
mainVect =c("a.","b.")
plot(Y[,1],Y[,2], pch=16, ylim=c(0,max(Y[,-1])), cex=0,xlab="",ylab="Density")
for(i in 2:4)
{
    points(Y[,1],Y[,i], pch=16, col=c("red","blue","orange")[i-1])
    # polygon(x=c(Ybar[,1],rev(Ybar[,1])), y=c(qgamma(p=0.05,shape=Ybar[,i],1),rev(qgamma(p=0.95,shape=Ybar[,i],1))), border=NA, col=grey(0.75,alpha=0.5))
    lines(Ybar[,1],Ybar[,i], col=c("red","blue","orange")[i-1])
}
legend("top",legend=c("X","Y","Z"),col=c("red","blue","orange"),lty=1,horiz=T,bty="n")
legend("topright",legend=mainVect[k],bty="n",cex=1.5)
k = k + 1
#
par(mfrow=c(1,1))

#
###
