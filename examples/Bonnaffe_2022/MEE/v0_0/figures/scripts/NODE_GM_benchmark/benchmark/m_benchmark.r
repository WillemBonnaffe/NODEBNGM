############
## main.R ##
############

## contact: Willem Bonnaffe (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

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

######################
## LOAD TIME SERIES ##
######################

## goal: load time series

## load data
TS = as.matrix(read.table("NODEBNGM/data/TS_3DLV.csv",sep=";",header=T))
TS = TS[20:50,]
TS[,1] = TS[,1]-min(TS[,1])
N = ncol(TS) - 1

## load ground truth
load("NODEBNGM/data/GT_3DLV.RData")
t_true = LV_GT$t_true
rhat_true = 1/LV_GT$Yhat_true*LV_GT$ddt.Yhat_true
ddx.rhat_true = LV_GT$ddx.rhat_true

## load results NODEBNGM
load("NODEBNGM/out/results_NODEBNGM.RData")

## load results NODE
load("NODE/out/results_NODE.RData")

## load results ODE2
load("ODE2/out/results_ODE2.RData")

## load results CCM
load("CCM/out/results_CCM.RData")

## list results
results = list("NODEBNGM"=results_NODEBNGM, "NODE"=results_NODE, "ODE2"=results_ODE2, "CCM"=results_CCM)

#
###

#############
## FIGURES ##
#############

##
## EFFECTS AND CONTRIBUTIONS

layout(mat = matrix(1:(N*3),nrow=3))
colVect = c("green","blue","red")
#
for(i in 1:N)
{
    ## per-capita growth rate
    plot(c(-100,100),y=c(0,0),type="l",lty=2,xlim=c(20,50)-1,ylim=c(-1,1),xlab="Time",ylab="Growth rate",cex.lab=1.5)
    for(k in 1:4) lines(results[[k]]$t,results[[k]]$Yhat_p[,i],col=colVect[i])
    lines(t_true,rhat_true[,i],col=colVect[i],lty=2)
    #
    ## effects
    plot(c(-100,100),y=c(0,0),type="l",lty=2,xlim=c(20,50)-1,ylim=c(-1,1)*3,cex=0,xlab="Time",ylab="Effects",cex.lab=1.5)
    #
    for(k in 1:4)
    {
        for(j in 1:N)
        {
            lines(results[[k]]$t,results[[k]]$ddx.Yhat_p[,j+(i-1)*N],col=colVect[j],lty=k)
        }
    }
    #
    for(j in 1:N) lines(t_true,ddx.rhat_true[,j+(i-1)*N],col=colVect[j],lty=2)
    #
    legend("bottom",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
    #
    ## geber
    plot(c(-100,100),y=c(0,0),type="l",lty=2,xlim=c(20,50)-1,ylim=c(-1,1)*1,cex=0,xlab="Time",ylab="Contributions",cex.lab=1.5)
    #
    for(k in 1:4)
    {
        for(j in 1:N)
        {
            lines(results[[k]]$t,results[[k]]$Geber_p[,j+(i-1)*N],col=colVect[j],lty=k)
        }
    }
    #
    legend("bottom",legend=colnames(TS)[-1],col=colVect,lty=1,horiz=T,bty="n")
}
#
par(mfrow=c(1,1))

##
## DIN

## compute mean Jacobians
J = list()
C = list()
MSq = function(x) mean(x^2,na.rm=T)
prop = function(x) x/sum(x)
for(k in 1:4) 
{
    J[[k]] = t(matrix(apply(results[[k]]$ddx.Yhat_p,2,mean,na.rm=T),ncol=N))
    C[[k]] = prop(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)))
    # C[[k]] = t(apply(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)),1,prop))
}

## visualise
par(mfrow=c(2,2))
#
for(k in 1:4) .plot.DIN(J[[k]],C[[k]],labels=colnames(TS)[-1])
#
par(mfrow=c(1,1))

#
###