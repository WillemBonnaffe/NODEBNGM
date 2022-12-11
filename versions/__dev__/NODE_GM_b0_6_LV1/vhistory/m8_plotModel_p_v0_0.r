#####################
## m7_visModel_p.r ##
#####################

## goal: visualise results of process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 08-04-2022 - created v0_0
## 15-04-2022 - updated dependencies

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## imports 
source("m1_loadData_o.r")
source("m2_loadModel_o.r")
source("m5_loadModel_p.r")

## load results
load(paste(pathToOut,"/","Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","Geber_p.RData",sep=""))

## output file
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))

#
###

########################
## PLOT PROCESS MODEL ##
########################

## goal: visualise predictions of process model

## figure predictions ##
col    = rainbow(N)
xlab   = c("","","Time")
ylab   = c("Dynamics","Effects","Contributions")
index  = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
legend = paste(colnames(TS)[-1]) 
par(mar=c(4,4,0,0),oma=c(1,1,1,1),cex.lab=1.25)
layout(cbind(c(1,2,3),c(4,5,6),c(7,8,9))[,1:min(3,N)])
for(i in 1:N)
{
	## attach predictions
	attach(    Yhat_p[[i]],warn.conflicts=F)
	attach(ddx.Yhat_p[[i]],warn.conflicts=F)
	attach(   Geber_p[[i]],warn.conflicts=F)
	#
	## training data
	source("m5_loadData_p.r")
	Y  = 1/Yhat_o[[i]][["E.Yhat_o"]]*ddt.Yhat_o[[i]][["E.ddt.Yhat_o"]] # Y_p[,i]
	#
	## dynamics
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(Y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=if(i==1)ylab[1]else"")
	points(nt,Y)
	polygon(c(nt,rev(nt)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col[i],alpha=0.2),border=NA)
	lines(nt,E.Yhat_p,col=adjustcolor(col[i],alpha=0.75))
	if(!is.null(index)) legend("topright",legend=index[1+(i-1)*(3)],bty="n",cex=1)
	#
	## effects
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=if(i==1)ylab[2]else"")
	for(j in 1:N) lines(nt,E.ddx.Yhat_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[2+(i-1)*(3)],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft" ,legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
	#
	## Geber 
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=if(i==1)ylab[3]else"")
	for(j in 1:N) lines(nt, E.Geber_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[3+(i-1)*(3)],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft" ,legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
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
