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
# TS = as.matrix(read.table("data/TS_3.csv",sep=";",header=T))
# TS[,1] = TS[,1]-min(TS[,1])
# N = ncol(TS) - 1
# TS[,-1][TS[,-1]<0.005] = 0.005

#
###

##############################
## FUNCTIONS NEURAL NETWORK ##
##############################

## activation functions
nn.pol = function(x) (x^2)
nn.relu = function(x) (x>0)*x
nn.sigmoid = function(x) 1/(1+exp(-x))
nn.exp = function(x) exp(x)

## nn.linear
## goal:
nn.linear = function(x,weights,I,O)
{
  weights = matrix(weights,nrow=I,ncol=O)
  return(x%*%weights)
}

## nn.fully_connected
nn.SLP = function(x,weights_1,weights_2,I,W,O)
{
  x = nn.linear(x,weights_1,I,W)
  x = nn.sigmoid(x)
  x = nn.linear(x,weights_2,W,O)
}

#
###

# ################
# ## DEBUG NODE ##
# ################
# 
# ##
# ## linear layer with relu activations
# O = 1
# I = 4
# input   = rnorm(I)
# weights = rnorm(I*O)
# print(nn.relu(nn.linear(input,weights,I,O)))
# 
# ##
# ## single layer perceptron
# I = 3
# W = 10
# O = 3
# input   = rnorm(I)
# weights_1 = rnorm(I*W)
# weights_2 = rnorm(W*O)
# print(nn.SLP(input,weights_1,weights_2,I,W,O))
# 
# ##
# ## polynomial linear model
# I = 3
# O = 1
# input    = rnorm(I)
# weights_1 = rnorm(I*O)
# weights_2 = rnorm((I^2)*O)
# print(nn.linear(input,weights_1,I,O) + nn.linear(as.vector(input%*%t(input)),weights_2,I^2,O))
# 
# #
# ###

################################
## FUNCTIONS DYNAMICAL MODELS ##
################################

## model.dYdt_AR
## goal:
n_Beta_AR = N^2
model.dYdt_AR = function(t, Y, Beta)
{
  dY = nn.linear(Y,Beta,N,N)
  return(as.vector(dY))
}

## model.dYdt_AR2
## goal:
n_Beta_AR2 = N^2 + N^3
model.dYdt_AR2 = function(t, Y, Beta)
{
  dY = nn.linear(Y,Beta[1:(N^2)],N,N) + nn.linear(as.vector(Y%*%t(Y)),Beta[-c(1:(N^2))],N^2,N)
  return(as.vector(dY))
}

## model.dYdt_ODE
n_Beta_ODE = N^2
model.dYdt_ODE = function(t, Y, Beta)
{
  dYdt = nn.linear(Y,Beta[1:(N^2)],N,N)*Y
  return(list(dYdt))
}

## model.dYdt_ODE2
n_Beta_ODE2 = N^2 + N^3
model.dYdt_ODE2 = function(t, Y, Beta)
{
  dYdt = (nn.linear(Y,Beta[1:(N^2)],N,N) + nn.linear(as.vector(Y%*%t(Y)),Beta[-c(1:(N^2))],N^2,N))*Y
  return(list(dYdt))
}

## model.dYdt_resNet
## goal:
I = N
W = 10
O = N
idx_weights_1 = 1:(I*W)
idx_weights_2 =   (I*W) + 1:(W*O)
idx_weights_3 =   (I*W) +   (W*O) + 1:(N*N)
n_Beta_resNet = I*W + W*O + N*N
model.dYdt_resNet = function(t, Y, Beta)
{
  x = nn.linear(Y,Beta[idx_weights_1],I,W)
  x = nn.exp(x)
  x = nn.linear(x,Beta[idx_weights_2],W,O)
  x = x + nn.linear(Y,Beta[idx_weights_3],N,N)
  x = x * Y
  return(as.vector(x))
}

## model.dYdt_NODE
## goal:
I = N
W = 10
O = N
idx_weights_1 = 1:(I*W)
idx_weights_2 =   (I*W) + 1:(W*O)
idx_weights_3 =   (I*W) +   (W*O) + 1:(N*N)
n_Beta_NODE   = I*W + W*O + N*N
model.dYdt_NODE = function(t, Y, Beta)
{
  x = nn.linear(Y,Beta[idx_weights_1],I,W)
  x = nn.exp(x)
  x = nn.linear(x,Beta[idx_weights_2],W,O)
  x = x + nn.linear(Y,Beta[idx_weights_3],N,N)
  x = x * Y
  return(list(as.vector(x)))
}

#
###

#############################
## FUNCTIONS MODEL FITTING ##
#############################

## model.predict_DE
## goal:
model.predict_DE = function(times, Y_0, Beta)
{
  Ybar = NULL
  Y = Y_0
  for(t in times)
  {
    Y = Y + model.dYdt(t,Y,Beta)
    Ybar = rbind(Ybar,Y)
  }
  Ybar = cbind(times,Ybar)
  return(Ybar)
}

## model.predict_ODE
## goal:
model.predict_ODE = function(times, Y_0, Beta)
{
  return(ode(y=Y_0, times=times, func=model.dYdt, parms=Beta, method="ode23"))
}

## model.dLogLik
## goal:
model.dLogLik = function(TS, Y_0, Beta, Sigma)
{
  TS_pred = model.predict(times=TS[,1], Y_0=Y_0, Beta=Beta)
  log_lik = sum(log(dnorm(log(TS[,-1]),log(TS_pred[,-1]),Sigma)))
  return(log_lik)
}

## model.dLogPrior
## goal:
model.dLogPrior = function(Y_0,Beta,Sigma)
{
  log_pri = sum(log(dunif(Y_0,0,10))) + 
    sum(log(dunif(Beta,-10,10))) + 
    sum(log(dnorm(log(Sigma),.5,.5)))
  return(log_pri)
}

## model.dLogPost
## goal:
model.dLogPost = function(TS, Y_0, Beta, Sigma)
{
  log_lik = model.dLogLik(TS, Y_0, Beta, Sigma)
  log_pri = model.dLogPrior(Y_0, Beta, Sigma)
  log_pos = log_lik + log_pri
  return(log_pos)
}

## model.dLogPost_wrapper
## goal:
model.dLogPost_wrapper = function(TS,Theta)
{
  Y_0   = exp(Theta[idx_Y_0])
  Beta  = Theta[idx_Beta]
  Sigma = exp(Theta[idx_Sigma])
  res   = model.dLogPost(TS, Y_0, Beta, Sigma)
  if(!is.nan(res) & !is.na(res))
  {
    return(res)
  }else
  {
    return(-Inf)
  }
}

## model.initiate
## goal:
model.initiate = function()
{
  check = F
  while(check == F)
  {
    Y_0   = TS[1,-1]
    Beta  = runif(n_Beta,-.001,.001)
    Sigma = runif(N,0,3)
    Theta_0 = c(log(Y_0),log(Sigma),Beta)
    if(model.dLogPost_wrapper(TS,Theta_0) > -Inf) 
    {
      check = T
    }
  }
  return(Theta_0)
}

## model.plot
## goal:
model.plot = function(TS,TS_pred)
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
    lines(TS_pred[,1],TS_pred[,i], col=colVect[i-1])
  }
  legend("top",legend=c("X","Y","Z"),col=colVect,lty=1,horiz=T,bty="n")
  legend("topright",legend=mainVect[k],bty="n",cex=1.5)
  k = k + 1
  #
  par(mfrow=c(1,1))
}

## model.predict
## goal:
## WIP - think about Ybar as predict function and model.predict as model.predict_wrapper
model.predict_wrapper = function(Theta)
{
  Y_0  = exp(Theta[idx_Y_0])
  Beta = Theta[idx_Beta]
  Ybar = model.predict(times=TS[,1], Y_0=Y_0, Beta=Beta)
  return(Ybar)
}

## model.dYdt
model.dYdt = model.dYdt_AR
model.predict = model.predict_DE
n_Beta = n_Beta_AR

## model properties
idx_Y_0   = 1:N
idx_Sigma =   N + 1:N
idx_Beta  =   N +   N + 1:n_Beta

## debug
Y_0   = as.numeric(TS[1,-1])
Beta  = runif(n_Beta,-.001,.001)
Sigma = rep(1,N)
TS_pred = model.predict(times = TS[,1], Y_0 = Y_0, Beta = Beta)
model.dLogPost(TS,Y_0,Beta,Sigma)
model.dLogPost_wrapper(TS,c(log(Y_0),log(Sigma),Beta))
model.dLogPost_wrapper(TS,model.initiate())

## benchmark
time = system.time(
  for(i in 1:1000)
  {
    model.dLogPost_wrapper(TS,model.initiate())
  }
)
print(time)

#
###

#################
## DEMC CHAINS ##
#################

## benchmark
timeVect = c(0,0)

## RCpp implementation of DEMCO
chainList = list()
Theta_0   = model.initiate()
tmax      = seq(10,15,1)
timeVect[1] = system.time(
  for(i in 1:length(tmax))
  {
    dTarget = function(x) model.dLogPost_wrapper(TS=TS[1:tmax[i],],x)
    chainList[[1]] = DEMCOpp(list(dTarget = dTarget,
                                 Theta_0 = Theta_0,
                                 gamma   = 2.38/sqrt(2*n_Beta),
                                 epsilon = 0.001,
                                 lambda  = 100,
                                 nIt     = 2000))[["chainList"]]
    Theta_0 = chainList.argmaxPost(chainList)
    model.plot(TS[1:tmax[i],],model.predict_wrapper(Theta_0))
    print(model.dLogPost_wrapper(TS=TS[1:tmax[i],], Theta=Theta_0 ))
  })[3]
MaP = chainList.argmaxPost(chainList)

## RCpp implementation of DEMC
chainList = list()
timeVect = c(0)
Theta_0 = MaP
timeVect[2] = system.time(
    for(i in 1:3)
    {
        # Theta_0 = initiate()
        # Theta_0 = Theta_0
        dTarget = function(x) model.dLogPost_wrapper(TS=TS,x)
        chainList[[i]] = DEMCpp(list(dTarget = dTarget,
                                     Theta_0 = Theta_0,
                                     gamma   = 2.38/sqrt(2*n_Beta),
                                     epsilon = 0.001,
                                     nIt     = 10000))[["chainList"]]
    })[3]
MaP = chainList.argmaxPost(chainList)
model.plot(TS,model.predict_wrapper(MaP))

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
Ybar_MaP = model.predict(times = TS[,1], Y_0 = Y_0, Beta = Beta)
#
## ensemble predictions
Ybar_ensemble = chainList.apply(chainList_thinned,f = function(x) model.predict(times = TS[,1], Y_0 = x[-1][idx_Y_0], Beta = x[-1][idx_Beta]))
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

