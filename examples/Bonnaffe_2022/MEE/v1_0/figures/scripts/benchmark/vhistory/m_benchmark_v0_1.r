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
.plot.DIN = function(effectsMat,weightsMat,labels,main="")
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
    plot(x=c(-1:1)*2,y=c(-1:1)*2,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main=main)
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
ddt.Yhat_true = LV_GT$ddt.Yhat_true
rhat_true = 1/LV_GT$Yhat_true*LV_GT$ddt.Yhat_true
ddx.rhat_true = LV_GT$ddx.rhat_true
Geber_true = cbind(ddt.Yhat_true,ddt.Yhat_true,ddt.Yhat_true) * ddx.rhat_true

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

## match time steps in ground truth to predictions
s = match(results[[1]]$t,t_true)
t_true = t_true[s]
ddt.Yhat_true = ddt.Yhat_true[s,] 
rhat_true = rhat_true[s,]
ddx.rhat_true = ddx.rhat_true[s,]
Geber_true = Geber_true[s,]

## remove extra time step in CCM
s = length(results$CCM$t)
results$CCM$t = results$CCM$t[-s]
results$CCM$Yhat_p = results$CCM$Yhat_p[-s,]
results$CCM$ddx.Yhat_p = results$CCM$ddx.Yhat_p[-s,]

#
###

##############
## ANALYSIS ##
##############

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

J_true = t(matrix(apply(ddx.rhat_true,2,mean,na.rm=T),ncol=N))
C_true = prop(t(matrix(apply(Geber_true,2,MSq),ncol=N)))

#
###

#############
## FIGURES ##
#############

##
## BENCHMARK

## graphical
pdf("out/benchmark.pdf",width=12,height=12)
#
## graphical parameters
layout(mat=matrix(1:(4*5),ncol=5))
color_vector = rainbow(4)
#
for(i in 1:N)
{
    ## per capita growth rate
    xlab = "Predicted"
    ylab = "Real"
    main = paste("P.c. growth rate of ",colnames(TS)[-1][i])
    plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(-1,1),ylim=c(-1,1),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
    #
    ## GT
    points(rhat_true[,i],rhat_true[,i],pch=16)
    #
    ## methods
    for(k in rev(1:4))
    {
        points(results[[k]]$Yhat_p[,i],rhat_true[,i],col=color_vector[k],pch=16)
    }
    #
    ## legend
    if (i==1) legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    
    ## effects
    for(j in 1:N)
    {
        y = x = ddx.rhat_true[,j+(i-1)*N]
        alpha = 0.5
        main = paste("Effect of ",colnames(TS)[-1][j]," on ", colnames(TS)[-1][i],sep="")
        xlab = "Predicted"
        ylab = "Real"
        plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(min(x)-alpha*(max(x)-min(x)),max(x)+alpha*(max(x)-min(x))),ylim=c(min(y)-alpha*(max(y)-min(y)),max(y)+alpha*(max(y)-min(y))),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
        #
        ## GT
        points(ddx.rhat_true[,j+(i-1)*N],ddx.rhat_true[,j+(i-1)*N],pch=16)
        #
        ## methods
        for(k in rev(1:4))
        {
        
            points(results[[k]]$ddx.Yhat_p[,j+(i-1)*N],ddx.rhat_true[,j+(i-1)*N],col=color_vector[k],pch=16)
        }
        #
        ## legend
        # legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    }
}
#
## dynamical interaction plot
for(k in 1:4) .plot.DIN(J[[k]],C[[k]],labels=colnames(TS)[-1],main=names(results)[k])
for(k in 1:4) .plot.DIN(J_true,C_true,labels=colnames(TS)[-1],main="Real")
#
par(mfrow=c(1,1))
#
dev.off()

##
## SUMMARY TABLE

## compute sum of squares
MSq = function(x) mean(x^2,na.rm=T)
res_growth = list()
res_effects = list()
res_contributions = list()
for(k in 1:4)
{
    res_growth[[k]]        = MSq(rhat_true - results[[k]]$Yhat_p)
    res_effects[[k]]       = MSq(ddx.rhat_true - results[[k]]$ddx.Yhat_p)
}
#
## barplot
par(mfrow=c(2,1))
barplot(unlist(res_growth),names.arg = names(results),ylab="MSq growth rate")
barplot(unlist(res_effects),names.arg = names(results),ylab="MSq effects")
par(mfrow=c(1,1))
#
## create table
benchmark_table = rbind(unlist(res_growth),unlist(res_effects))

#
###