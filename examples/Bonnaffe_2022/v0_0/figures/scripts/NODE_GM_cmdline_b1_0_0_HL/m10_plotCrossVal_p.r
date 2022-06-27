#########################
## m8_plotCrossVal_p.r ##
#########################

## goal: visualise results of cross validation

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## imports 
source("m5_loadData_p.r")
source("m6_loadModel_p.r")

## load results
load(paste(pathToOut,"/","crossVal_p.RData",sep=""))

## output file
pdf(paste(pathToOut,"/fig_crossVal_p.pdf",sep=""))

#
###

##########
## PLOT ##
##########

## goal: visualise cross validation 

## figure cross validation ##
index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
colVect = rainbow(2,start=0.6,end=0.9)
# par(mfrow=c(min(3,N),1),mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
layout(cbind(c(1,2),c(3,4),c(5,6))[,1:min(3,N)])
for(i in 1:N)
{
	## unpack
	crossVal = crossVal_p[[i]]
	#
	## plot
	x    = crossVal[,"sd"]
	# y    = crossVal[,c("logLik_l","logLik_t")]
	# sd.y = crossVal[,c("sd.logLik_l","sd.logLik_t")]
    xlab = "Std. prior"
	ylab = if(i == 1) "Log likelihood" else ""
	#
	## training line
	col  = colVect[1]
	y    = crossVal[,"logLik_l"] 
	sd.y = crossVal[,"sd.logLik_l"]
	plot(x,rep(0,length(x)),ylim=c(min(y-sd.y)-0.5*abs(min(y-sd.y)),max(y+sd.y)+0.5*abs(max(y+sd.y))),type="l",lty=3,xlab="",ylab=ylab)
	lines(x,y,col=col,lwd=2)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
    ## legend
	legend("bottom",legend=c("Training","Testing")[1],col=colVect[1],lty=1,bty="n",cex=1.5,horiz=T,lwd=2)
	legend("topright",legend=index[1+(i-1)*2],bty="n",cex=1.5)
    #
	## testing line
	col  = colVect[2]	
	y    = crossVal[,"logLik_t"]
	sd.y = crossVal[,"sd.logLik_t"]
	plot(x,rep(0,length(x)),ylim=c(min(y-sd.y)-0.5*abs(min(y-sd.y)),max(y+sd.y)+0.5*abs(max(y+sd.y))),type="l",lty=3,xlab=xlab,ylab=ylab)
	lines(x,y,col=col,lwd=2)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
	## legend
	legend("bottom",legend=c("Training","Testing")[2],col=colVect[2],lty=1,bty="n",cex=1.5,horiz=T,lwd=2)
	legend("topright",legend=index[2+(i-1)*2],bty="n",cex=1.5)
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
