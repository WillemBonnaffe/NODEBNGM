############
## main.R ##
############

##############
## INITIATE ##
##############

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

## load ground truth
load("data/GT_3DLV.RData")
ddx.rhat_true_t = LV_GT$t_true
ddx.rhat_true = LV_GT$ddx.rhat_true

#
###

#####################
## FUNCTIONS UTILS ##
#####################

utils.pairwise = function(x) as.vector(x%*%t(x))
utils.ddx.pairwise = function(x) apply(diag(length(x)),2,function(y) as.vector(y%*%t(x) + x%*%t(y)))

#
###

##############################
## FUNCTIONS NEURAL NETWORK ##
##############################

## activation functions
nn.pol = function(x) (x^2)
nn.ddx.pol = function(x) 2*x
nn.relu = function(x) (x>0)*x
nn.ddx.relu = function(x) (x>0)
nn.sigmoid = function(x) 1/(1+exp(-x))
nn.ddx.sigmoid = function(x) nn.sigmoid(x)*(1-nn.sigmoid(x))
nn.exp = function(x) exp(x)
nn.ddx.exp = function(x) exp(x)

## nn.linear
## goal:
nn.linear = function(x,weights,I,O)
{
  weights = matrix(weights,nrow=I,ncol=O)
  return(x%*%weights)
}

## nn.ddx.linear
## goal:
nn.ddx.linear = function(x,weights,I,O)
{
  weights = matrix(weights,nrow=I,ncol=O)
  # return(rep(1,I)%*%weights)
  return(weights)
}

## nn.fully_connected
nn.SLP = function(x0,weights_1,weights_2,I,W,O)
{
  x1 = nn.linear(x0,weights_1,I,W)
  x2 = nn.sigmoid(x1)
  x3 = nn.linear(x2,weights_2,W,O)
  return(x3)
}
nn.ddx.SLP = function(x0,weights_1,weights_2,I,W,O)
{
  ## forward pass
  x1 = nn.linear(x0,weights_1,I,W)
  x2 = nn.sigmoid(x1)
  x3 = nn.linear(x2,weights_2,W,O)

  ## derivatives
  dx1dx0 = nn.ddx.linear(x0,weights_1)
  dx2dx1 = nn.ddx.sigmoid(x1)
  dx3dx2 = nn.ddx.linear(x2,weights_2)
  return(dx1dx0*dx2dx1*dx3*dx2)
}

#
###

################################
## FUNCTIONS DYNAMICAL MODELS ##
################################

## model.dYdt_AR
## goal:
n_Beta_AR = function(I,O) I*O
model.dYdt_AR = function(t, Y, Beta, I, O)
{
  dYdt = nn.linear(Y,Beta,I,O)*Y
  return(as.vector(dYdt))
}
model.ddx.dYdt_AR = function(t, Y, Beta, I, O)
{
  ddx.dYdt = nn.ddx.linear(Y,Beta,I,O)
  return(ddx.dYdt)
}

## model.dYdt_AR2
## goal:
n_Beta_AR2 = function(I,O) I*O + I*I*O
model.dYdt_AR2 = function(t, Y, Beta, I, O)
{
  weights_1 = Beta[ c(1:(I*O))]
  weights_2 = Beta[-c(1:(I*O))]
  #
  x0 = Y
  x1 = nn.linear(x0,weights_1,I,O)
  x2 = utils.pairwise(x0)
  x3 = nn.linear(x2,weights_2,I*I,O)*Y
  x4 = x1 + x3
  x5 = x4*Y
  #
  return(as.vector(x5))
}
model.ddx.dYdt_AR2 = function(t, Y, Beta, I, O)
{
  weights_1 = Beta[ c(1:(I*O))]
  weights_2 = Beta[-c(1:(I*O))]
  #
  x0 = Y
  x1 = nn.linear(x0,weights_1,I,O)
  x2 = utils.pairwise(x0)
  x3 = nn.linear(x2,weights_2,I*I,O)*Y
  x4 = x1 + x3
  x5 = x4*Y
  #
  dx1dx0 = nn.ddx.linear(x0,weights_1,I,O)
  dx2dx0 = utils.ddx.pairwise(x0)
  dx3dx2 = nn.ddx.linear(x2,weights_2,I*I,O)
  dx4dx0 = dx1dx0 + t(dx2dx0)%*%dx3dx2
  #
  # x = Y 
  # w_1 = Beta[1:(I*O)]
  # w_2 = Beta[-c(1:(I*O))]
  # ddx.dYdt = nn.ddx.linear(x,w_1,I,O) + 
  #   utils.ddx.pairwise(x) %*% nn.ddx.linear(utils.pairwise(x),w_2,I*I,O)
  #
  return(dx4dx0)
}

## debug
I = 3
O = 4
Y = rnorm(I)
weights_1 = rnorm(I*O)
weights_2 = rnorm(I*I*O)
## WIP - to ground truth


## model.dYdt_ODE
n_Beta_ODE = function(I,O) I*O
model.dYdt_ODE = function(t, Y, Beta, I, O)
{
  dYdt = nn.linear(Y,Beta,I,O)*Y
  return(list(dYdt))
}

## model.dYdt_ODE2
n_Beta_ODE2 = function(I,O) I*I + I*I*O
model.dYdt_ODE2 = function(t, Y, Beta, I, O)
{
  dYdt = (nn.linear(Y,Beta[1:(I*O)],I,O) + nn.linear(as.vector(Y%*%t(Y)),Beta[-c(1:(I*O))],I*I,O))*Y
  return(list(dYdt))
}

## model.dYdt_resNet
## goal:
n_Beta_resNet = function(I,W,O) I*W + W*O + O*O
model.dYdt_resNet = function(t, Y, Beta, I, W=10, O)
{
  idx_weights_1 = 1:(I*W)
  idx_weights_2 =   (I*W) + 1:(W*O)
  idx_weights_3 =   (I*W) +   (W*O) + 1:(O*O)
  x = nn.linear(Y,Beta[idx_weights_1],I,W)
  x = nn.exp(x)
  x = nn.linear(x,Beta[idx_weights_2],W,O)
  x = x + nn.linear(Y,Beta[idx_weights_3],O,O)
  x = x * Y
  return(as.vector(x))
}

## model.dYdt_NODE
## goal:
n_Beta_NODE = function(I,W,O) I*W + W*O + O*O
model.dYdt_NODE = function(t, Y, Beta, I, W=10, O)
{
  idx_weights_1 = 1:(I*W)
  idx_weights_2 =   (I*W) + 1:(W*O)
  idx_weights_3 =   (I*W) +   (W*O) + 1:(O*O)
  x = nn.linear(Y,Beta[idx_weights_1],I,W)
  x = nn.exp(x)
  x = nn.linear(x,Beta[idx_weights_2],W,O)
  x = x + nn.linear(Y,Beta[idx_weights_3],O,O)
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
  return(ode(y=Y_0, times=times, func=model.dYdt, parms=Beta, method="ode45"))
}

## model.predict
## goal:
## WIP - think about Ybar as predict function and model.predict as model.predict_wrapper
model.predict_wrapper = function(Theta)
{
  Y_0  = exp(Theta[idx_Y_0()])
  Beta = Theta[idx_Beta()]
  Ybar = model.predict(times=TS[,1], Y_0=Y_0, Beta=Beta)
  return(Ybar)
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
  Y_0   = exp(Theta[idx_Y_0()])
  Sigma = exp(Theta[idx_Sigma()])
  Beta  = Theta[idx_Beta()]
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
  # mainVect =c("a.","b.")
  colVect = rainbow(N)
  plot(TS[,1],TS[,2], pch=16, ylim=c(0,max(TS[,-1])), cex=0,xlab="",ylab="Density")
  for(i in 2:(N+1))
  {
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    lines(TS_pred[,1],TS_pred[,i], col=colVect[i-1])
  }
  legend("top",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
  # legend("topright",legend=mainVect[k],bty="n",cex=1.5)
  # k = k + 1
  #
  par(mfrow=c(1,1))
}

#
###

####################
## INITIATE MODEL ##
####################

## model properties
N = ncol(TS) - 1

## define model.dYdt
model.dYdt = function(t, Y, Beta) model.dYdt_AR2(t, Y, Beta, I=N, O=N)
model.ddx.dYdt = function(t, Y, Beta) model.ddx.dYdt_AR2(t, Y, Beta, I=N, O=N)
n_Beta = n_Beta_AR2(I=N, O=N)
model.predict = model.predict_DE
output = "out_AR2"

## fixed input
idx_Y_0   = function() 1:N
idx_Sigma = function()  N + 1:N
idx_Beta  = function()  N +   N + 1:n_Beta
n_params = length(idx_Y_0()) + length(idx_Sigma()) + length(idx_Beta())

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
tmax      = seq(10,30,1)
timeVect[1] = system.time(
  for(i in 1:length(tmax))
  {
    dTarget = function(x) model.dLogPost_wrapper(TS=TS[1:tmax[i],],x)
    chainList[[1]] = DEMCOpp(list(dTarget = dTarget,
                                 Theta_0 = Theta_0,
                                 gamma   = 2.38/sqrt(2*n_params),
                                 epsilon = 0.001,
                                 lambda  = 100,
                                 nIt     = 1000))[["chainList"]]
    Theta_0 = chainList.argmaxPost(chainList)
    model.plot(TS[1:tmax[i],],model.predict_wrapper(Theta_0))
    print(model.dLogPost_wrapper(TS=TS[1:tmax[i],], Theta=Theta_0 ))
  })[3]
MaP = chainList.argmaxPost(chainList)

## RCpp implementation of DEMC
chainList = list()
timeVect = c(0)
Theta_0 = MaP
n_chains = 3
timeVect[2] = system.time(
    for(i in 1:n_chains)
    {
        # Theta_0 = initiate()
        # Theta_0 = Theta_0
        dTarget = function(x) model.dLogPost_wrapper(TS=TS,x)
        chainList[[i]] = DEMCpp(list(dTarget = dTarget,
                                     Theta_0 = Theta_0,
                                     gamma   = 2.38/sqrt(2*n_params),
                                     epsilon = 0.001,
                                     nIt     = 10000))[["chainList"]]
    })[3]
MaP = chainList.argmaxPost(chainList)
model.plot(TS,model.predict_wrapper(MaP))

## store chains
system(paste("rm -r ",output,sep=""))
system(paste("mkdir ",output,sep=""))
chainList.write(chainList,output)

#
###

######################
## DEMC DIAGNOSTICS ##
######################

## load chains
chainList = chainList.read(output)
# timeVect = c(312.362,125.786,102.705)

## untransform
for(l in 1:length(chainList))
{
  chainList[[l]][,-1][,idx_Y_0()]  = exp(chainList[[l]][,-1][,idx_Y_0()])
  chainList[[l]][,-1][,idx_Sigma()] = exp(chainList[[l]][,-1][,idx_Sigma()])
  chainList[[l]][,-1][,idx_Beta()]  = chainList[[l]][,-1][,idx_Beta()]
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

##
## FIT

png(paste(output,"/fit.png",sep=""))
#
## MaP prediction
MaP  = as.numeric(chainList.argmaxPost(chainList_thinned))
Y_0  = MaP[idx_Y_0()]
Beta = MaP[idx_Beta()]
Ybar_MaP = model.predict(times = TS[,1], Y_0 = Y_0, Beta = Beta)
#
## ensemble predictions
Ybar_ensemble = chainList.apply(chainList_thinned,f = function(x) model.predict(times = TS[,1], Y_0 = x[-1][idx_Y_0()], Beta = x[-1][idx_Beta()]))
Ybar_mean = matrix(Ybar_ensemble[[1]],ncol=N+1)
Ybar_q05  = matrix(Ybar_ensemble[[2]],ncol=N+1)
Ybar_q95  = matrix(Ybar_ensemble[[3]],ncol=N+1)
#
## plot
par(mfrow=c(1,1),mar=c(5,5,0,0),oma=c(0,0,1,1))
#
# k = 1
# mainVect =c("a.","b.")
colVect = rainbow(N)
plot(TS[,1],TS[,2], pch=16, ylim=c(0,max(TS[,-1])), cex=0,xlab="",ylab="Density")
for(i in 2:4)
{
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    polygon(x=c(Ybar_q05[,1],rev(Ybar_q05[,1])), y=c(Ybar_q05[,i],rev(Ybar_q95[,i])), border=NA, col=grey(0.75,alpha=0.5))
    lines(Ybar_MaP[,1],Ybar_MaP[,i], col=colVect[i-1])
}
legend("top",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
# legend("topright",legend=mainVect[k],bty="n",cex=1.5)
# k = k + 1
#
par(mfrow=c(1,1))
#
dev.out()

##
## Effects

##
ddx.rhat_pred = t(apply(Ybar_mean, 1, function(x) model.ddx.dYdt(t = x[1], Y = x[-1], Beta = MaP[idx_Beta()])))

##
par(mfrow=c(2,N))
#
for(i in 1:N)
{
  plot(1:10,xlim=c(20,50)-1,ylim=c(-1,1)*3.5,cex=0)
  #
  t = TS[,1] + 19
  lines(t,ddx.rhat_pred[,1]*0,lty=2)
  lines(t,ddx.rhat_pred[,1+(i-1)*N],col="green")
  lines(t,ddx.rhat_pred[,2+(i-1)*N],col="blue")
  lines(t,ddx.rhat_pred[,3+(i-1)*N],col="red")
  #
  t = ddx.rhat_true_t
  lines(t,ddx.rhat_true[,1+(i-1)*N],col="green",lty=2)
  lines(t,ddx.rhat_true[,2+(i-1)*N],col="blue",lty=2)
  lines(t,ddx.rhat_true[,3+(i-1)*N],col="red",lty=2)
}
par(mfrow=c(1,1))

#
###

