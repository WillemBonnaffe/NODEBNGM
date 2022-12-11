############
## main.R ##
############

## contact: Willem Bonnaffe (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

## libraries
library(deSolve)

## source functions
source("functions/f_HBM_v1_0.r")

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
t_true = LV_GT$t_true
rhat_true = 1/LV_GT$Yhat_true*LV_GT$ddt.Yhat_true
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

#####################
## FUNCTIONS PLOTS ##
#####################

## .plot.DIN
## goal: plot the dynamical interaction network of the system
# effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
# weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
# labels     - vector - vector of the names of the variables in the matrix
.plot.DIN = function(effectsMat,weightsMat,labels)
{
    ## dimensions
    N = dim(effectsMat)[1]
    
    ## scale effects and contributions
    effectsMat = (effectsMat>0)*1
    # weightsMat = weightsMat/sum(weightsMat) # proportion of total change
    
    ## angles
    theta = seq(0,2*pi,(2*pi)/N)
    x = cos(theta)
    y = sin(theta)
    x_ = cos(theta+(2*pi*0.05))
    y_ = sin(theta+(2*pi*0.05))
    x__ = 1.25*cos(theta+(2*pi*0.025))
    y__ = 1.25*sin(theta+(2*pi*0.025))
    
    ## plot interactions
    plot(x=c(-1:1)*2,y=c(-1:1)*2,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
    for(i in 1:N)
    {
        points(x[i],y[i],pch=16)
        for(j in 1:N)
        {
            color_ = if(effectsMat[i,j]>0){"green"}else{"red"}
            # points(x__[i],y__[i],cex=30/N)
            text(x__[i],y__[i],labels=labels[i])
            if(weightsMat[i,j]*10>0)
            {
                arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
            }
        }
    }
}

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
  return(weights)
}

#
###

################################
## FUNCTIONS DYNAMICAL MODELS ##
################################

## model.dYdt_ODE
n_Beta_ODE = function(I,O) I*O
model.dYdt_ODE = function(t, Y, Beta, I, O)
{
  dYdt = nn.linear(Y,Beta,I,O)*Y
  return(list(dYdt))
}
model.ddx.dYdt_ODE = function(t, Y, Beta, I, O)
{
  ddx.dYdt = nn.ddx.linear(Y,Beta,I,O)
  return(ddx.dYdt)
}

## model.dYdt_ODE2
n_Beta_ODE2 = function(I,O) I*I + I*I*O
model.dYdt_ODE2 = function(t, Y, Beta, I, O)
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
  return(list(x5))
}
model.ddx.dYdt_ODE2 = function(t, Y, Beta, I, O)
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
  return(dx4dx0)
}

## model.dYdt_NODE
## goal:
n_Beta_NODE = function(I,W,O) I*W + W*O + I*O
model.dYdt_NODE = function(t, Y, Beta, I, W=10, O)
{
  weights_1 = Beta[1:(I*W)]
  weights_2 = Beta[  (I*W) + 1:(W*O)]
  weights_3 = Beta[  (I*W) +   (W*O) + 1:(I*O)]
  #
  x0 = Y
  x1 = nn.linear(x0,weights_1,I,W)
  x2 = nn.exp(x1)
  x3 = nn.linear(x2,weights_2,W,O)
  x4 = nn.linear(x0,weights_3,I,O)
  x5 = x3 + x4
  x6 = x5 * x0
  #
  return(list(x6))
}
model.ddx.dYdt_NODE = function(t, Y, Beta, I, W=10, O)
{
  ## debug >>
  # I = 3
  # W = 10
  # O = 5
  # Y = rnorm(I)
  # Beta = rnorm(n_Beta_NODE(I,W,O))
  ## debug <<
  weights_1 = Beta[1:(I*W)]
  weights_2 = Beta[  (I*W) + 1:(W*O)]
  weights_3 = Beta[  (I*W) +   (W*O) + 1:(I*O)]
  #
  x0 = Y
  x1 = nn.linear(x0,weights_1,I,W)
  x2 = nn.exp(x1)
  x3 = nn.linear(x2,weights_2,W,O)
  x4 = nn.linear(x0,weights_3,I,O)
  x5 = x3 + x4
  x6 = x5 * x0
  #
  dx1dx0 = nn.ddx.linear(x0,weights_1,I,W)
  dx2dx1 = nn.ddx.exp(x1)
  dx3dx2 = nn.ddx.linear(x2,weights_2,W,O)
  dx4dx0 = nn.ddx.linear(x0,weights_3,I,O)
  dx5dx0 = t(apply(dx1dx0,1,function(x)x*dx2dx1))%*%dx3dx2 + dx4dx0
  #
  return(dx5dx0)
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
    Y = Y + unlist(model.dYdt(t,Y,Beta))
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
  colVect = rainbow(N)
  plot(TS[,1],TS[,2], pch=16, ylim=c(0,max(TS[,-1])), cex=0,xlab="",ylab="Density")
  for(i in 2:(N+1))
  {
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    lines(TS_pred[,1],TS_pred[,i], col=colVect[i-1])
  }
  legend("top",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
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
model.dYdt = function(t, Y, Beta) model.dYdt_ODE2(t, Y, Beta, I=N, O=N)
model.ddx.dYdt = function(t, Y, Beta) model.ddx.dYdt_ODE2(t, Y, Beta, I=N, O=N)
n_Beta = n_Beta_ODE2(I=N, O=N)
model.predict = model.predict_ODE
output = "out"

## fixed input
idx_Y_0   = function() 1:N
idx_Sigma = function()  N + 1:N
idx_Beta  = function()  N +   N + 1:n_Beta
n_params = length(idx_Y_0()) + length(idx_Sigma()) + length(idx_Beta())

## check
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

## RCpp implementation of DEMCO
t_train      = 1:20
t_val        = 21:30
n_iterations = 60
n_chains     = 30
chainList    = list()
timeVect     = rep(0,n_chains)
for(l in 1:n_chains)
{
  timeVect[l] = system.time({

      ## initiate
      dTarget = function(x) - model.dLogPost_wrapper(TS=TS[c(t_train),],x)
      Theta_0 = Theta  = model.initiate()
      dTarget_0 = - dTarget(Theta_0)

      ## pre-training
      try(expr = {
        res     = optim(Theta,dTarget,method="BFGS",control = list("trace"=1,"REPORT"=1,"maxit"=n_iterations))
        Theta = res$par
      })
      
      # ## initiate
      # dTarget = function(x) - model.dLogPost_wrapper(TS=TS[c(t_train,t_val),],x)
      # 
      # ## training
      # try(expr = {
      #   res     = optim(Theta,dTarget,method="BFGS",control = list("trace"=1,"REPORT"=1,"maxit"=n_iterations))
      #   Theta = res$par
      # })

      ## store in change
      Theta_f = Theta
      dTarget_f = - dTarget(Theta_f)
      chainList[[l]] = rbind(c(dTarget_0,Theta_0),c(dTarget_f,Theta_f))

      # ## visualise
      # model.plot(TS,model.predict_wrapper(Theta))
      # print(model.dLogPost_wrapper(TS=TS[t_train,], Theta=Theta))
      # print(model.dLogPost_wrapper(TS=TS[t_val,], Theta=Theta))

    })[3]
}

## store chains
system(paste("rm -r ",output,sep=""))
system(paste("mkdir ",output,sep=""))
chainList.write(chainList,output)
write(x = timeVect, file = paste(output,"/runtimes.txt",sep=""))

#
###

######################
## CHAIN PROCESSING ##
######################

## load chains
chainList = chainList.read(output)
chainList_thinned = chainList

## untransform
for(l in 1:length(chainList_thinned))
{
  chainList_thinned[[l]][,-1][,idx_Y_0()]  = exp(chainList_thinned[[l]][,-1][,idx_Y_0()])
  chainList_thinned[[l]][,-1][,idx_Sigma()] = exp(chainList_thinned[[l]][,-1][,idx_Sigma()])
  chainList_thinned[[l]][,-1][,idx_Beta()]  = chainList_thinned[[l]][,-1][,idx_Beta()]
}

## group chains
chainList_thinned = chainList.unlist(chainList_thinned)
chainList_thinned = chainList_thinned[!is.na(chainList_thinned[,1]),]
chainList_thinned = chainList_thinned[chainList_thinned[,1]>quantile(x=chainList_thinned[,1], p=0.75),]
chainList_thinned = list(chainList_thinned)

#
###

#############
## FIGURES ##
#############

##
## FIT

## ensemble predictions
Ybar_ensemble = chainList.apply(chainList_thinned,f = function(x) model.predict(times = TS[,1], Y_0 = x[-1][idx_Y_0()], Beta = x[-1][idx_Beta()]))
Ybar_mean = matrix(Ybar_ensemble[[1]],ncol=N+1)
Ybar_q05  = matrix(Ybar_ensemble[[2]],ncol=N+1)
Ybar_q95  = matrix(Ybar_ensemble[[3]],ncol=N+1)

## avoid na on test set
s = which(is.na(Ybar_ensemble$f_ensemble[,ncol(Ybar_ensemble$f_ensemble)]))
chainList_thinned[[1]] = chainList_thinned[[1]][-s,]

## MaP prediction
MaP  = as.numeric(chainList.argmaxPost(chainList_thinned))
Y_0  = MaP[idx_Y_0()]
Beta = MaP[idx_Beta()]
Ybar_MaP = model.predict(times = TS[,1], Y_0 = Y_0, Beta = Beta)

## figure
pdf(paste(output,"/fit.pdf",sep=""))
#
par(mfrow=c(1,1),mar=c(5,5,0,0),oma=c(0,0,1,1))
#
colVect = c("green","blue","red")
plot(TS[,1],TS[,2], pch=16, ylim=c(0,3), cex=0,xlab="Time",ylab="Density")
for(i in 2:4)
{
    points(TS[,1],TS[,i], pch=16, col=colVect[i-1])
    polygon(x=c(Ybar_q05[,1],rev(Ybar_q05[,1])), y=c(Ybar_q05[,i],rev(Ybar_q95[,i])), border=NA, col=grey(0.75,alpha=0.5))
    lines(Ybar_MaP[,1],Ybar_MaP[,i], col=colVect[i-1],lty=2)
    lines(Ybar_mean[,1],Ybar_mean[,i], col=colVect[i-1])
}
legend("top",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
#
par(mfrow=c(1,1))
#
dev.off()

##
## EFFECTS

## compute per-capita growth rate effects and contributions
rhat_pred     = t(apply(Ybar_MaP, 1, function(x) 1/x[-1] * unlist(model.dYdt(t = x[1], Y = x[-1], Beta = MaP[idx_Beta()]))))
ddx.rhat_pred = t(apply(Ybar_MaP, 1, function(x) model.ddx.dYdt(t = x[1], Y = x[-1], Beta = MaP[idx_Beta()])))
geber_pred    = t(apply(Ybar_MaP, 1, function(x) model.ddx.dYdt(t = x[1], Y = x[-1], Beta = MaP[idx_Beta()]) * unlist(model.dYdt(t = x[1], Y = x[-1], Beta = MaP[idx_Beta()]))))

## figure
pdf(paste(output,"/effects.pdf",sep=""))
#
layout(mat = matrix(1:(N*3),nrow=3))
colVect = c("green","blue","red")
#
for(i in 1:N)
{
    ## per-capita growth rate
    plot(1:10,xlim=c(20,50)-1,ylim=c(-1,1),cex=0,xlab="Time",ylab="Growth rate",cex.lab=1.5)
    #
    t = TS[,1] + 19
    lines(t,rhat_pred[,i],col=colVect[i])
    #
    t = t_true
    lines(t,rhat_true[,i],col=colVect[i],lty=2)
    #
    ## effects
    plot(1:10,xlim=c(20,50)-1,ylim=c(-1,1)*3,cex=0,xlab="Time",ylab="Effects",cex.lab=1.5)
    #
    t = TS[,1] + 19
    lines(t,ddx.rhat_pred[,1]*0,lty=2)
    lines(t,ddx.rhat_pred[,1+(i-1)*N],col="green")
    lines(t,ddx.rhat_pred[,2+(i-1)*N],col="blue")
    lines(t,ddx.rhat_pred[,3+(i-1)*N],col="red")
    #
    t = t_true
    lines(t,ddx.rhat_true[,1+(i-1)*N],col="green",lty=2)
    lines(t,ddx.rhat_true[,2+(i-1)*N],col="blue",lty=2)
    lines(t,ddx.rhat_true[,3+(i-1)*N],col="red",lty=2)
    #
    legend("bottom",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
    #
    ## geber
    plot(1:10,xlim=c(20,50)-1,ylim=c(-1,1),cex=0,xlab="Time",ylab="Contributions",cex.lab=1.5)
    #
    t = TS[,1] + 19
    lines(t,geber_pred[,1]*0,lty=2)
    lines(t,geber_pred[,1+(i-1)*N],col="green")
    lines(t,geber_pred[,2+(i-1)*N],col="blue")
    lines(t,geber_pred[,3+(i-1)*N],col="red")
    #
    legend("bottom",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
}
par(mfrow=c(1,1))
#
dev.off()

#
###

############
## OUTPUT ##
############

## goal: store results

## save predictions
results_ODE2 = list("t"=TS[,1]+19, "Yhat_p"=rhat_pred, "ddx.Yhat_p"=ddx.rhat_pred, "Geber_p"=geber_pred)
save(results_ODE2,file=paste(output,"/results_ODE2.RData",sep=""))

#
###