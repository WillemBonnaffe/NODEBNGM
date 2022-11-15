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
TS = as.matrix(read.table("data/TS_3DLV.csv",sep=";",header=T))
TS = TS[20:50,]
TS[,1] = TS[,1]-min(TS[,1])
N = ncol(TS) - 1

# ## load data
# TS = as.matrix(read.table("data/TS_1.csv",sep=";",header=T))
# TS[,1] = TS[,1]-min(TS[,1])
# N = ncol(TS) - 1
# TS[,-1][TS[,-1]<0.005] = 0.005

#
###

##################
## INITIATE ODE ##
##################

## model properties
idx_Y_0   = 1:N
idx_Sigma =   N + 1:N
idx_Beta  =   N +   N + 1:(N^2)

## dynamics
dYdt = function(t, Y, Beta)
{
  J = matrix(Beta,ncol=length(Y)) 
  dYdt = (J%*%Y)*Y
  return(list(dYdt))
}

## predictive function
Ybar = function(times, Y_0, dYdt, Beta)
{
  return(ode(y=Y_0, times=times, func=dYdt, parms=Beta, method="ode23"))
}

## posterior
dLogPost = function(TS, Ybar, Y_0, Beta, Sigma)
{
  ## global parameters
  N = ncol(TS) - 1
  
  ## predictive function
  TS_pred = Ybar(times=TS[,1], Y_0=Y_0, dYdt=dYdt, Beta=Beta)
  
  ## calculate posterior
  log_lik = sum(log(dnorm(log(TS[,-1]),log(TS_pred[,-1]),Sigma)))
  log_pri = sum(log(dunif(Beta,-10,10))) + sum(log(dnorm(log(Sigma),.5,.5)))
  log_pos = log_lik + log_pri
  
  # ## debug
  # print(log_pos)
  # print(log_lik)
  # print(log_pri)
  
  ## terminate
  return(log_pos)
}

## wrapper
dLogPost_wrapper = function(TS,paramVect)
{
  Y_0   = exp(paramVect[idx_Y_0])
  Sigma = exp(paramVect[idx_Sigma])
  Beta  = paramVect[idx_Beta]
  res   = dLogPost(TS, Ybar, Y_0, Beta, Sigma)
  if(!is.nan(res) & !is.na(res)){return(res)}else{return(-Inf)}
}

## DEMC wrapper
dTarget = function(x)
{
  return(dLogPost_wrapper(TS=TS,x))
}

## initiate function
initiate = function()
{
  check = F
  while(check == F)
  {
    Y_0   = TS[1,-1]
    Beta  = runif(N^2,-.001,.001)
    Sigma = runif(N,0,3)
    Theta_0 = c(log(Y_0),log(Sigma),Beta)
    # names(Beta_0) = c("V1","V2","V3","V4")
    if(dTarget(Theta_0) > -Inf){check = T}
  }
  return(Theta_0)
}

## plot function
model.plot = function(TS,Ybar_MaP)
{
  par(mfrow=c(1,1),mar=c(5,5,0,0),oma=c(0,0,1,1))
  #
  k = 1
  mainVect =c("a.","b.")
  colVect = rainbow(N)
  plot(TS[,1],TS[,2], pch=16, ylim=c(0,max(TS[,-1])), cex=0,xlab="",ylab="Density")
  for(i in 2:(N+1))
  {
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    # polygon(x=c(Ybar_q05[,1],rev(Ybar_q05[,1])), y=c(Ybar_q05[,i],rev(Ybar_q95[,i])), border=NA, col=grey(0.75,alpha=0.5))
    lines(Ybar_MaP[,1],Ybar_MaP[,i], col=colVect[i-1])
  }
  legend("top",legend=c("X","Y","Z"),col=colVect,lty=1,horiz=T,bty="n")
  legend("topright",legend=mainVect[k],bty="n",cex=1.5)
  k = k + 1
  #
  par(mfrow=c(1,1))
}

model.predict = function(MaP)
{
  Y_0  = exp(MaP[idx_Y_0])
  Beta = MaP[idx_Beta]
  Ybar_MaP = Ybar(times = TS[,1], Y_0 = Y_0, dYdt = dYdt, Beta = Beta)
  return(Ybar_MaP)
}

## debug
Y_0   = TS[1,-1]
Beta  = runif(N^2,-.001,.001)
Sigma = rep(1,N)
TS_pred = Ybar(times = TS[,1], Y_0 = Y_0, dYdt = dYdt, Beta = Beta)
dLogPost(TS,Ybar,Y_0,Beta,Sigma)
dTarget(c(log(Y_0),log(Sigma),Beta)) # checks out
dTarget(initiate())

#
###

################
## ODE CHAINS ##
################

## benchmark
timeVect = c(0,0)

## RCpp implementation of DEMCO
chainList = list()
Theta_0 = initiate()
tmax    = seq(10,30,2)
# tmax = c(10,20,30,40,50,60)
timeVect[1] = system.time(
  for(i in 1:length(tmax))
  {
    dTarget = function(x) dLogPost_wrapper(TS=TS[1:tmax[i],],x)
    chainList[[1]] = DEMCOpp(list(dTarget = dTarget,
                                 Theta_0 = Theta_0,
                                 gamma   = 2.38/sqrt(2*N),
                                 epsilon = 0.001,
                                 lambda  = 100,
                                 nIt     = 1000))[["chainList"]]
    Theta_0 = chainList.argmaxPost(chainList)
    model.plot(TS[1:tmax[i],],model.predict(Theta_0))
  })[3]
MaP = chainList.argmaxPost(chainList)

## RCpp implementation of DEMC
chainList = list()
timeVect = c(0)
Theta_0 = MaP
timeVect[2] = system.time(
    for(i in 1:1)
    {
        # Theta_0 = initiate()
        # Theta_0 = Theta_0
        dTarget = function(x) dLogPost_wrapper(TS=TS,x)
        chainList[[i]] = DEMCpp(list(dTarget = dTarget,
                                     Theta_0 = Theta_0,
                                     gamma   = 2.38/sqrt(2*N),
                                     epsilon = 0.001,
                                     nIt     = 10000))[["chainList"]]
    })[3]
MaP = chainList.argmaxPost(chainList)
model.plot(TS,model.predict(MaP))

## store chains
chainList.write(chainList,"out")

#
###

######################
## DEMC DIAGNOSTICS ##
######################

## load chains
chainList = chainList.read("out")
# timeVect = c(312.362,125.786,102.705)

## untransform
for(l in 1:length(chainList))
{
  chainList[[l]][,-1][,idx_Y_0]  = exp(chainList[[l]][,-1][,idx_Y_0])
  chainList[[l]][,-1][,idx_Sigma] = exp(chainList[[l]][,-1][,idx_Sigma])
  chainList[[l]][,-1][,idx_Beta]  = chainList[[l]][,-1][,idx_Beta]
}

## traces MC
burnin = c(1:5000)
chainList.tracePlot(chainList)
chainList.tracePlot(chainList.thin(chainList.burn(chainList,burnin)))

## plot posterior distribution DEMC
chainList.postPlot(chainList.thin(chainList.burn(chainList,burnin)),1000)

## ac plots MC
par(mfrow=c(3,4))
chainList.acPlot(chainList.burn(chainList,burnin))
par(mfrow=c(1,1))

## summary table
summaryTab = chainList.summaryTab(chainList.thin(chainList.burn(chainList,burnin)),returnSignif=F)[["estimates"]][-1,]
ESSTab = chainList.ESS(chainList.thin(chainList.burn(chainList,burnin)))
#
# write.table(comparisonTab,"out_ODE/comparisonTab.txt",sep=" & ",quote=F, row.names=F)

## thin
chainList_thinned = chainList.thin(chainList.burn(chainList,burnin))

#
###

#############
## FIGURES ##
#############

png("out/fit.png")
#
## MaP prediction
MaP  = as.numeric(chainList.argmaxPost(chainList_thinned))
Y_0  = MaP[idx_Y_0]
Beta = MaP[idx_Beta]
Ybar_MaP = Ybar(times = TS[,1], Y_0 = Y_0, dYdt = dYdt, Beta = Beta)
#
## ensemble predictions
Ybar_ensemble = chainList.apply(chainList_thinned,f = function(x) Ybar(times = TS[,1], Y_0 = x[-1][idx_Y_0], dYdt = dYdt, Beta = x[-1][idx_Beta]))
Ybar_q05 = matrix(Ybar_ensemble[[2]],ncol=N+1)
Ybar_q95 = matrix(Ybar_ensemble[[3]],ncol=N+1)
#
## plot
par(mfrow=c(1,1),mar=c(5,5,0,0),oma=c(0,0,1,1))
#
k = 1
mainVect =c("a.","b.")
colVect = rainbow(N)
plot(TS[,1],TS[,2], pch=16, ylim=c(0,max(TS[,-1])), cex=0,xlab="",ylab="Density")
for(i in 2:4)
{
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    polygon(x=c(Ybar_q05[,1],rev(Ybar_q05[,1])), y=c(Ybar_q05[,i],rev(Ybar_q95[,i])), border=NA, col=grey(0.75,alpha=0.5))
    lines(Ybar_MaP[,1],Ybar_MaP[,i], col=colVect[i-1])
}
legend("top",legend=c("X","Y","Z"),col=colVect,lty=1,horiz=T,bty="n")
legend("topright",legend=mainVect[k],bty="n",cex=1.5)
k = k + 1
#
par(mfrow=c(1,1))
#
dev.out()


##
J = matrix(Beta,ncol=N)
image(J)

#
###

