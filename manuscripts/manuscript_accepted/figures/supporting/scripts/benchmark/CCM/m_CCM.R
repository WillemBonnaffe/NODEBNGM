#############
## m_CCM.r ## 
#############

## goal: perform convergent cross mapping analysis of artificial tri-trophic system

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

###############
## FUNCTIONS ##
###############

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

##############
## INITIATE ##
##############

## libraries
library(rEDM)

#
###

######################
## LOAD TIME SERIES ##
######################

## load data
TS = as.matrix(read.table("data/TS_3DLV.csv",sep=";",header=T))
TS = TS[20:50,]
# TS[,-1] = log(TS[,-1])# [20:50,]
TS[,1] = TS[,1]-min(TS[,1])
N = ncol(TS)-1

## format for smaps
TS = data.frame(TS)
lib = c(1, 20)
pred = c(1,31)
cols = c("G", "B", "R")
targets = c("G", "B", "R")

## load ground truth
load("data/GT_3DLV.RData")
t_true = LV_GT$t_true
Yhat_true = LV_GT$Yhat_true
rhat_true = 1/LV_GT$Yhat_true*LV_GT$ddt.Yhat_true
ddx.rhat_true = LV_GT$ddx.rhat_true

## output
output = "out"

#
###

##########
## MAIN ##
##########

## goal: explain the dynamics of each variables as a function of other variables

predictions_list = list()
effects_list = list()
contributions_list = list()
runtimes = rep(0,length(targets))
for(i in 1:length(targets))
{
    
    ## update 
    print(paste(i,"/",length(targets),sep=""))
    
    ## select target
    target = targets[i]
    
    ## S-map
    runtimes[i] = system.time({
        block_smap_output = block_lnlp(TS, 
                                        lib = lib, 
                                        pred = pred, 
                                        columns = cols, 
                                        target_column = target, 
                                        method = "s-map", 
                                        theta = 2, 
                                        stats_only = FALSE, 
                                        first_column_time = TRUE, 
                                        save_smap_coefficients = TRUE, 
                                        silent = TRUE)
    })[3]
        
    ## extract S-map coefs
    smap_coeffs = block_smap_output$smap_coefficients[[1]]
    
    ## make prediction
    predictions = block_smap_output$model_output[[1]]
    
    ## approximate dynamics
    dxdt = apply(rbind(predictions,predictions[nrow(predictions),]),2,diff)
    
    ## compute effects
    effects = smap_coeffs[-c(1:2)]
    
    ## compute contributions
    contributions = effects * dxdt
    
    ## store
    predictions_list[[i]] = predictions
    effects_list[[i]] = effects
    contributions_list[[i]] = contributions
    
}

#
###

##############
## ANALYSIS ##
##############

## compute mean Jacobian matrix
J = NULL
C = NULL
MSq = function(x) mean(x^2,na.rm=T)
prop = function(x) x/sum(x)
for(i in 1:length(targets))
{
    
    ## extract effects and contributions
    effects = effects_list[[i]]
    contributions = contributions_list[[i]]
    
    ## compute mean effects and contributions
    J_ = apply(effects,2,mean,na.rm=T)
    C_ = prop(apply(contributions,2,MSq))
    
    ## build Jacobian and contributions
    J = rbind(J,J_)
    C = rbind(C,C_) 
}

## compute per-capita growth rate
growth_rate_list = list()
for(i in 1:length(targets))
{
    t = predictions_list[[i]]$t
    predictions = predictions_list[[i]]$Predictions
    rhat_pred = 1/predictions * diff(c(predictions,NA))
    growth_rate_list[[i]] = data.frame(cbind("t"=t,"growth_rate"=rhat_pred))
}

#
###

#############
## FIGURES ##
#############

## goal: visualise predictions and effects obtained from S-maps

##
## FIT

## figure
pdf(paste(output,"/fit.pdf",sep=""))
#
## graphical parameters
color_vector = c("green","blue","red")
#
plot(1:10,xlim=c(20,50)-1,ylim=c(0,3),cex=0,xlab="Time",ylab="Density",cex.lab=1.5)
#
for (i in 1:length(targets))
{
    ## extract predictions and smap_coeffs
    predictions = predictions_list[[i]]
    growth_rate = growth_rate_list[[i]]
    effects = effects_list[[i]]
    contributions = contributions_list[[i]]
    #
    ## predictions
    points(predictions$t+19, predictions$Observations, col=color_vector[i],pch=16)
    lines(predictions$t+19, predictions$Predictions, col=color_vector[i])
    lines(t_true, Yhat_true[,i], lty=2, col=color_vector[i])
    #
    legend("top",legend=colnames(TS)[-1],col=color_vector,lty=1,horiz=T,bty="n")
 
}
#
par(mfrow=c(1,1))
#
dev.off()

##
## EFFECTS

pdf(paste(output,"/effects.pdf",sep=""))
#
## graphical parameters
par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
layout(mat=matrix(1:9,ncol=3))
color_vector = c("green","blue","red")
index  = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
#
for (i in 1:length(targets))
{
    
    ## extract predictions and smap_coeffs
    predictions = predictions_list[[i]]
    growth_rate = growth_rate_list[[i]]
    effects = effects_list[[i]]
    contributions = contributions_list[[i]]
    
    ## growth rate
    plot(growth_rate$t, growth_rate$growth_rate, type="l", ylim=c(-1,1), col=color_vector[i], xlab="",ylab=if(i==1)"Growth rate"else"",cex.lab=1.5)
    lines(t_true-19, rhat_true[,i], lty=2, col=color_vector[i])
    if(!is.null(index)) legend("topright",legend=index[1+(i-1)*(3)],bty="n",cex=1.5)
    
    ## effects
    plot(c(-1000,1000), c(0,0), xlim=c(0,max(predictions$t)), ylim=c(-1,1)*3, type="l", lty=2, xlab="", ylab=if(i==1)"Effects"else"",cex.lab=1.5)
    for(j in 1:length(cols))
    {
        lines(predictions$t, effects[, j],   col = color_vector[j])
        lines(t_true-19, ddx.rhat_true[,j+(i-1)*N],   col = color_vector[j], lty=2)
    }
    legend("bottom",legend=colnames(TS)[-1],col=color_vector,lty=1,horiz=T,bty="n")
    if(!is.null(index)) legend("topright",legend=index[2+(i-1)*(3)],bty="n",cex=1.5)
    
    ## contributions
    plot(c(-1000,1000), c(0,0), xlim=c(0,max(predictions$t)), ylim=c(-1,1)*1, type="l", lty=2, xlab="Time", ylab=if(i==1)"Contributions"else"",cex.lab=1.5)
    for(j in 1:length(cols))
    {
        lines(predictions$t, contributions[, j],   col = color_vector[j])
    }
    #
    legend("bottom",legend=colnames(TS)[-1],col=color_vector,lty=1,horiz=T,bty="n")
    if(!is.null(index)) legend("topright",legend=index[3+(i-1)*(3)],bty="n",cex=1.5)
}
#
par(mfrow=c(1,1))
#
dev.off()

##
## DYNAMICAL INTERACTION NETWORK
pdf(paste(output,"/DIN.pdf",sep=""))
.plot.DIN(J,C,colnames(TS)[-1])
dev.off()

#
###

############
## OUTPUT ##
############

## format results
t = predictions_list[[1]]$t + 19
Yhat_p = matrix(unlist(growth_rate_list),ncol=2*N)[,seq(2,N*2,2)]
ddx.Yhat_p = matrix(unlist(effects_list),ncol=N*N)
Geber_p = matrix(unlist(contributions_list),ncol=N*N)

## save results
results_CCM = list("t"=t,"Yhat_p"=Yhat_p,"ddx.Yhat_p"=ddx.Yhat_p,"Geber_p"=Geber_p)
save(results_CCM,file=paste(output,"/results_CCM.RData",sep=""))

## save runtimes
write(runtimes,paste(output,"/runtimes.txt",sep=""))

#
###
