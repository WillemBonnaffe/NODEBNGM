#######################
## m8_plotModel_p.r ##
#######################

## goal: visualise results of process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0
## 15-04-2022 - updated dependencies
## 25-04-2022 - created v0_1
##            - updated function repository
## 18-05-2022 - created v0_2
##            - embellished figures
## 26-05-2022 - created v0_3
##            - renamed m8_plotModel_p.r

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## imports 
source("m1_loadData_o.r")
source("m2_loadModel_o.r")
source("m5_loadData_p.r")
source("m6_loadModel_p.r")

## load results
load(paste(pathToOut,"/","Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","Geber_p.RData",sep=""))

## load ground truth
load("data/GT_3DLV.RData")
attach(LV_GT)

## output file
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))

#
###

########################
## PLOT PROCESS MODEL ##
########################

## goal: visualise predictions of process model

## figure predictions ##
col    = rev(rainbow(N,start=0.6,end=0.9))
xlab   = c("","","Time")
ylab   = c("P.c. growth rate","Effects","Contributions")
index  = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
legend = paste(colnames(TS)[-1]) 
par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
layout(cbind(c(1,2,3),c(4,5,6),c(7,8,9))[,1:min(3,N)])
for(i in 1:N)
{
	## predictions
	E.Yhat_p       = apply(Yhat_p[[i]],2,mean)
	q05.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.05)
	q95.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.95)
	E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p[[i]],2,mean),ncol=nrow(X)))
	q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
	q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
	E.Geber_p      = t(matrix(apply(Geber_p[[i]],2,mean),ncol=nrow(X)))
	q05.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
	q95.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
	#
	## training data
	Y  = E.Y_p[,i]
	#
	## dynamics
    x = nt
    y = Y
	plot(x,rep(0,length(x)),ylim=c(-1,1)*max(abs(y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=if(i==1)ylab[1]else"")
    # axis(1,at=seq(min(x),max(x),(max(x)-min(x))/5),labels=round(seq(min(x),max(x),(max(x)-min(x))/5)),lwd=0)
    # axis(2,at=seq(min(y),max(y),(max(y)-min(y))/5),labels=round(seq(min(y),max(y),(max(y)-min(y))/5),2),lwd=0)
	points(x,y,pch=16,col=adjustcolor("black",0.75))
	polygon(c(nt,rev(nt)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col[i],alpha=0.2),border=NA)
	lines(nt,E.Yhat_p,col=adjustcolor(col[i],alpha=0.75),lwd=2)
	if(!is.null(index)) legend("topright",legend=index[1+(i-1)*(3)],bty="n",cex=1.5)
    lines(t_true,1/Yhat_true[,i]*ddt.Yhat_true[,i],lty=2,col=adjustcolor("black",0.75),lwd=2)
	#
	## effects
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=if(i==1)ylab[2]else"")
	for(j in 1:N) lines(nt,E.ddx.Yhat_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
	for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
	if(!is.null(index))  legend("topright",legend=index[2+(i-1)*(3)],bty="n",cex=1.5)
	if(!is.null(legend)) legend("bottom" ,legend=legend,lty=1,col=adjustcolor(col,alpha=0.75),bty="n",horiz=T,lwd=2)
    for(j in 1:N) lines(t_true,ddx.rhat_true[,(1:3) + (i-1)*3][,j],lty=2,col=adjustcolor(col[j],alpha=0.75),lwd=2)
	#
	## Geber 
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=if(i==1)ylab[3]else"")
	for(j in 1:N) lines(nt, E.Geber_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
	for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
	if(!is.null(index))  legend("topright",legend=index[3+(i-1)*(3)],bty="n",cex=1.5)
	if(!is.null(legend)) legend("bottom" ,legend=legend,lty=1,col=adjustcolor(col,alpha=0.75),bty="n",horiz=T,lwd=2)
}

#
###

###############
## TERMINATE ##
###############

## terminate
dev.off()

#
###

############
## OUTPUT ##
############

## format results
Yhat_p_ = NULL
ddx.Yhat_p_ = NULL
Geber_p_ = NULL
for(i in 1:N)
{
    ## predictions
    E.Yhat_p       = apply(Yhat_p[[i]],2,mean)
    E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p[[i]],2,mean),ncol=nrow(X)))
    E.Geber_p      = t(matrix(apply(Geber_p[[i]],2,mean),ncol=nrow(X)))
    
    ## store 
    Yhat_p_ = cbind(Yhat_p_,E.Yhat_p)
    ddx.Yhat_p_ = cbind(ddx.Yhat_p_,E.ddx.Yhat_p)
    Geber_p_ = cbind(Geber_p_,E.Geber_p)
}

## save results
results_NODEBNGM = list("t"=nt,"Yhat_p"=Yhat_p_,"ddx.Yhat_p"=ddx.Yhat_p_,"Geber_p"=Geber_p_)
save(results_NODEBNGM,file=paste(pathToOut,"/","results_NODEBNGM.RData",sep=""))
