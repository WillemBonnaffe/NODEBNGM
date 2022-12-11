############
## main.R ##
############

## contact: Willem Bonnaffe (w.bonnaffe@gmail.com)

## to do:
## perform cross validation on train set NODEBNGM
## run optimisation CCM NODE and ODE2 only on train set
## update runtimes for all methods

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

## load results
load("NODEBNGM/out/results_NODEBNGM.RData")
load("NODE/out/results_NODE.RData")
load("ODE2/out/results_ODE2.RData")
load("CCM/out/results_CCM.RData")

## remove extra time step in CCM
s = length(results_CCM$t)
results_CCM$t = results_CCM$t[-s]
results_CCM$Yhat_p = results_CCM$Yhat_p[-s,]
results_CCM$ddx.Yhat_p = results_CCM$ddx.Yhat_p[-s,]
results_CCM$Geber_p = results_CCM$Geber_p[-s,]

## load ground truth
load("NODEBNGM/data/GT_3DLV.RData")
t_true = LV_GT$t_true
ddt.Yhat_true = LV_GT$ddt.Yhat_true
rhat_true = 1/LV_GT$Yhat_true*LV_GT$ddt.Yhat_true
ddx.rhat_true = LV_GT$ddx.rhat_true
Geber_true = cbind(ddt.Yhat_true,ddt.Yhat_true,ddt.Yhat_true) * ddx.rhat_true

## match time steps in ground truth to predictions
s = match(results_NODEBNGM$t,t_true)
t_true = t_true[s]
ddt.Yhat_true = ddt.Yhat_true[s,] 
rhat_true = rhat_true[s,]
ddx.rhat_true = ddx.rhat_true[s,]
Geber_true = Geber_true[s,]

## list results
GT = list("t"=t_true, "Yhat_p"=rhat_true, "ddx.Yhat_p"=ddx.rhat_true, "Geber_p"=Geber_true)
results = list("GT"=GT,"NODEBNGM"=results_NODEBNGM, "NODE"=results_NODE, "ODE2"=results_ODE2, "CCM"=results_CCM)

## split training and test set
s = 1:20
results_train = results_test = results
for(k in 1:5)
{
    ## train set
    results_train[[k]]$t = results[[k]]$t[s]
    results_train[[k]]$Yhat_p = results[[k]]$Yhat_p[s,]
    results_train[[k]]$ddx.Yhat_p = results[[k]]$ddx.Yhat_p[s,]
    
    ## test set
    results_test[[k]]$t = results[[k]]$t[-s]
    results_test[[k]]$Yhat_p = results[[k]]$Yhat_p[-s,]
    results_test[[k]]$ddx.Yhat_p = results[[k]]$ddx.Yhat_p[-s,]
}

## load runtimes
runtimes = list()
runtimes$GT = 0
runtimes$NODEBNGM = mean(unlist(as.vector(read.table("NODEBNGM/out/runtimes_o.txt"))))*N + mean(unlist(as.vector(read.table("NODEBNGM/out/runtimes_p.txt"))))*N
runtimes$NODE = mean(unlist(as.vector(read.table("NODE/out/runtimes.txt"))))
runtimes$ODE2 = mean(unlist(as.vector(read.table("ODE2/out/runtimes.txt"))))
runtimes$CCM = mean(as.numeric(read.table("CCM/out/runtimes.txt")))*N

#
###

#############
## FIGURES ##
#############


##
## BENCHMARK FIGURE TRAIN SET

## select train or test set
results = results_train

## compute mean Jacobians
J = list()
C = list()
MSq = function(x) mean(x^2,na.rm=T)
prop = function(x) x/sum(x)
for(k in 1:5) 
{
    J[[k]] = t(matrix(apply(results[[k]]$ddx.Yhat_p,2,mean,na.rm=T),ncol=N))
    C[[k]] = prop(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)))
    # C[[k]] = t(apply(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)),1,prop))
}

## graphical
pdf("out/benchmark_train.pdf",width=12,height=12)
#
## graphical parameters
layout(mat=matrix(1:(4*3),ncol=3))
color_vector = adjustcolor(c(adjustcolor("black",alpha=0),rainbow(4)),alpha=0.75)
#
for(i in 1:N)
{
    ## per capita growth rate
    xlab = "Ground truth"
    ylab = "Predicted"
    main = paste("P.c. growth rate of ",colnames(TS)[-1][i])
    plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(-1,1),ylim=c(-1,1),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
    #
    ## methods
    for(k in rev(1:5))
    {
        points(results[[1]]$Yhat_p[,i],results[[k]]$Yhat_p[,i],col=color_vector[k],pch=16)
    }
    #
    ## legend
    if (i==1) legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    
    ## effects
    for(j in 1:N)
    {
        y = x = results$GT$ddx.Yhat_p[,j+(i-1)*N]
        alpha = 0.5
        main = paste("Effect of ",colnames(TS)[-1][j]," on ", colnames(TS)[-1][i],sep="")
        xlab = "Ground truth"
        ylab = "Predicted"
        plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(min(x)-alpha*(max(x)-min(x)),max(x)+alpha*(max(x)-min(x))),ylim=c(min(y)-alpha*(max(y)-min(y)),max(y)+alpha*(max(y)-min(y))),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
        #
        ## methods
        for(k in rev(1:5))
        {
        
            points(results[[1]]$ddx.Yhat_p[,j+(i-1)*N],results[[k]]$ddx.Yhat_p[,j+(i-1)*N],col=color_vector[k],pch=16)
        }
        #
        ## legend
        # legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    }
}
#
# ## dynamical interaction plot
# for(k in 2:5) .plot.DIN(J[[k]],C[[k]],labels=colnames(TS)[-1],main=names(results)[k])
# for(k in 2:5) .plot.DIN(J[[1]],C[[1]],labels=colnames(TS)[-1],main="Ground truth")
#
par(mfrow=c(1,1))
#
dev.off()

##
## BENCHMARK FIGURE TEST SET

## select train or test set
results = results_test

## compute mean Jacobians
J = list()
C = list()
MSq = function(x) mean(x^2,na.rm=T)
prop = function(x) x/sum(x)
for(k in 1:5) 
{
    J[[k]] = t(matrix(apply(results[[k]]$ddx.Yhat_p,2,mean,na.rm=T),ncol=N))
    C[[k]] = prop(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)))
    # C[[k]] = t(apply(t(matrix(apply(results[[k]]$Geber_p,2,MSq),ncol=N)),1,prop))
}

## graphical
pdf("out/benchmark_test.pdf",width=12,height=12)
#
## graphical parameters
layout(mat=matrix(1:(4*3),ncol=3))
color_vector = adjustcolor(c(adjustcolor("black",alpha=0),rainbow(4)),alpha=0.75)
#
for(i in 1:N)
{
    ## per capita growth rate
    xlab = "Ground truth"
    ylab = "Predicted"
    main = paste("P.c. growth rate of ",colnames(TS)[-1][i])
    plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(-1,1),ylim=c(-1,1),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
    #
    ## methods
    for(k in rev(1:5))
    {
        points(results[[1]]$Yhat_p[,i],results[[k]]$Yhat_p[,i],col=color_vector[k],pch=16)
    }
    #
    ## legend
    if (i==1) legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    
    ## effects
    for(j in 1:N)
    {
        y = x = results$GT$ddx.Yhat_p[,j+(i-1)*N]
        alpha = 0.5
        main = paste("Effect of ",colnames(TS)[-1][j]," on ", colnames(TS)[-1][i],sep="")
        xlab = "Ground truth"
        ylab = "Predicted"
        plot(c(-1,1)*10,y=c(-1,1)*10,type="l",lty=2,xlim=c(min(x)-alpha*(max(x)-min(x)),max(x)+alpha*(max(x)-min(x))),ylim=c(min(y)-alpha*(max(y)-min(y)),max(y)+alpha*(max(y)-min(y))),xlab=xlab,ylab=ylab,cex.lab=1.5,main=main)
        #
        ## methods
        for(k in rev(1:5))
        {
            
            points(results[[1]]$ddx.Yhat_p[,j+(i-1)*N],results[[k]]$ddx.Yhat_p[,j+(i-1)*N],col=color_vector[k],pch=16)
        }
        #
        ## legend
        # legend("topleft", legend=names(results), col=color_vector, pch=16, bty="n")
    }
}
#
## dynamical interaction plot
# for(k in 2:5) .plot.DIN(J[[k]],C[[k]],labels=colnames(TS)[-1],main=names(results)[k])
# for(k in 2:5) .plot.DIN(J[[1]],C[[1]],labels=colnames(TS)[-1],main="Ground truth")
#
par(mfrow=c(1,1))
#
dev.off()

##
## SUMMARY TABLE

## compute sum of squares
MSq = function(x) mean(x^2,na.rm=T)
r2  = function(y,yhat) 1 - mean((y-yhat)^2)/mean((mean(y)-y)^2)
res_growth_train = list()
res_effects_train = list()
res_growth_test = list()
res_effects_test = list()
for(k in 1:5)
{
    # ## r-squared
    # res_growth_train[[k]]  = r2(results_train[[1]]$Yhat_p,     results_train[[k]]$Yhat_p)
    # res_effects_train[[k]] = r2(results_train[[1]]$ddx.Yhat_p, results_train[[k]]$ddx.Yhat_p)
    # res_growth_test[[k]]  = r2(results_test[[1]]$Yhat_p,       results_test[[k]]$Yhat_p)
    # res_effects_test[[k]] = r2(results_test[[1]]$ddx.Yhat_p,   results_test[[k]]$ddx.Yhat_p)
    #
    ## mean squared error
    res_growth_train[[k]]  = MSq(results_train[[1]]$Yhat_p - results_train[[k]]$Yhat_p)
    res_effects_train[[k]] = MSq(results_train[[1]]$ddx.Yhat_p - results_train[[k]]$ddx.Yhat_p)
    res_growth_test[[k]]  = MSq(results_test[[1]]$Yhat_p - results_test[[k]]$Yhat_p)
    res_effects_test[[k]] = MSq(results_test[[1]]$ddx.Yhat_p - results_test[[k]]$ddx.Yhat_p)
}
#
## barplot
pdf("out/benchmark_overview.pdf",width=12,height=12)
#
layout(mat=rbind(1:5,rep(6,5),rep(7,5),rep(8,5)))
color_vector = adjustcolor(c("blue","red"),alpha=0.5)
#
for(k in 1:5) .plot.DIN(J[[k]],C[[k]],labels=colnames(TS)[-1],main=names(results)[k])
#
par(mar=rep(5,4))
#
y = rbind(unlist(res_growth_train),unlist(res_growth_test))
barplot(y,names.arg=names(results),ylab="MSE growth rate",col=color_vector,beside=T,border=NA,cex.lab=2,cex.names=1.5)
legend("topleft",legend=c("train","test"),col=color_vector,pch=15,bty="n",cex=2)
#
y = rbind(unlist(res_effects_train),unlist(res_effects_test))
barplot(y,names.arg=names(results),ylab="MSE effects",col=color_vector,beside=T,border=NA,cex.lab=2,cex.names=1.5)
legend("topleft",legend=c("train","test"),col=color_vector,pch=15,bty="n",cex=2)
#
y = unlist(runtimes)
barplot(y,names.arg=names(results),ylab="Runtimes (s)",border=NA,cex.lab=2,cex.names=1.5)
#
par(mfrow=c(1,1))
#
dev.off()
#
## create table
# benchmark_table = rbind(unlist(res_growth_train),unlist(res_growth_test),unlist(res_effects_train),unlist(res_effects_test))

#
###