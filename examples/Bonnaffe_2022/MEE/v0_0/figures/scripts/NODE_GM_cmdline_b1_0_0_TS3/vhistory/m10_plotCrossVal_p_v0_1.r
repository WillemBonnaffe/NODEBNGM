#########################
## m8_plotCrossVal_p.r ##
#########################

## goal: visualise results of cross validation

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 15-04-2022 - created v0_0 from m7_visModel_p.r
##            - specific plot for cross validation
## 26-05-2022 - created v0_1
##            - adapted plot to work with marginal likelihood

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
par(mfrow=c(min(3,N),1),mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
for(i in 1:N)
{
	## unpack
	crossVal = crossVal_p[[i]]
	#
	## plot
	x    = crossVal[,"w"]
	y    = crossVal[,c("logMarLik_l","logMarLik_t")]
	sd.y = crossVal[,c("sd.logMarLik_l","sd.logMarLik_t")]
	xlab = if(i == N) "Weight" else ""
	ylab = "Log Marginal likelihood"
	plot(x,rep(0,length(x)),ylim=c(min(y-sd.y)-0.5*abs(min(y-sd.y)),max(y+sd.y)+0.5*abs(max(y+sd.y))),type="l",lty=3,xlab=xlab,ylab=ylab)
	#
	## training line
	col  = colVect[1]
	y    = crossVal[,"logMarLik_l"] 
	sd.y = crossVal[,"sd.logMarLik_l"]
	lines(x,y,col=col)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
	## testing line
	col  = colVect[2]	
	y    = crossVal[,"logMarLik_t"]
	sd.y = crossVal[,"sd.logMarLik_t"]
	lines(x,y,col=col)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
	## legend
	legend("bottom",legend=c("Training","Testing"),col=colVect,lty=1,bty="n",cex=1.5,horiz=T)
	legend("topright",legend=index[i],bty="n",cex=1.5)
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
