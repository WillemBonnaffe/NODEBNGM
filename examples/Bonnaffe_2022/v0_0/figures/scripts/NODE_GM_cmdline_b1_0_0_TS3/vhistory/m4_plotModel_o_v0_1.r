######################
## m4_plotModel_o.r ##
######################

## goal: visualise results of observation model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0
## 15-04-2022 - created v0_1
##            - simplified code
##            - renamed to plotModel_o

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load data 
source("m1_loadData_o.r")

## load results
load(paste(pathToOut,"/","Yhat_o.RData",sep=""))
load(paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))

## output file
pdf(paste(pathToOut,"/fig_predictions_o.pdf",sep=""))

#
###

#################################
## VISUALISE OBSERVATION MODEL ##
#################################

## goal: visualise observation model 

## figure
col   = rainbow(N)
xlab  = c("","Time")
ylab  = c("Y(t)","dY/dt")
index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
par(mar=c(4,4,0,0),oma=c(1,1,1,1),cex.lab=1.25)
layout(cbind(c(1,2),c(3,4),c(5,6))[,1:min(3,N)])
for (i in 1:N)
{
	## predictions 
    E.Yhat_o       = apply(Yhat_o[[i]],2,mean)
    q05.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.05)
    q95.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.95)
    E.ddt.Yhat_o   = apply(ddt.Yhat_o[[i]],2,mean)
    q05.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.05)
    q95.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.95)
    #
	## training data
	t = X_o 
	Y = Y_o[,i] 
	nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)
	#
	## visualise interpolated response
	plot(t,rep(0,length(t)),ylim=c(min(Y),max(Y))+0.2*c(-min(Y),max(Y)),type="l",lty=3,xlab=xlab[1],ylab=if(i == 1)ylab[1]else"")
	polygon(c(nt,rev(nt)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
	lines(nt,E.Yhat_o,col=adjustcolor(col[i],0.75))
	points(t,Y) 
	if(!is.null(index)) legend("topright",legend=index[1+(i-1)*2],bty="n",cex=1)
	#
	## visualise temporal derivative
	plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*1.5,type="l",lty=3,xlab=xlab[2],ylab=if(i == 1)ylab[2]else"")
	polygon(c(nt,rev(nt)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
	lines(nt,E.ddt.Yhat_o,col=adjustcolor(col[i],0.75))
	if(!is.null(index)) legend("topright",legend=index[2+(i-1)*2],bty="n",cex=1)
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
