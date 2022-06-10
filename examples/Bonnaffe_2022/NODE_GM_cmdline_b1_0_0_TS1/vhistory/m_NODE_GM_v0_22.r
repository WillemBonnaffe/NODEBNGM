############
## main.r ##
############

## goal: load time series and apply the NODERMT analysis 

## update:
## 10-03-2021 - created version 0.0
## 11-03-2021 - created version 0.1
##            - implemented figures for the paper 
## 29-03-2021 - created version 0.2
##            - combined figures for the paper
## 25-05-2021 - created version 0.3
##            - modularised code to work with any single time series file
## 18-11-2021 - crated version 0_5
##            - cleaned code
## 22-11-2021 - cleaned code
## 22-11-2021 - created v0_6
##            - removed unused sections
## 22-11-2021 - crated v0_7
##            - transferred main functions to function file
## 23-11-2021 - created v0_8
##            - implemented scaling of the process model
## 26-11-2021 - created v0_9
##            - implemented a fitting method that does not require anchor sampling
## 29-11-2021 - created v0_10
##            - experimented further with code
## 30-11-2021 - created v0_11
##            - implemented Geber method
##            - created v0_12
##            - implemented main function
## 17-02-2022 - created v0_13
##            - cleaned code
##            - improved visualisations
## 18-02-2022 - created v0_14
##            - changed fitting function
##            - code summary function (do not calculate quantities on the fly)
##            - separate functions for fitting the observation and process model
## 23-02-2022 - created v0_15
##            - moved sampling and plotting code to function folder (f_NODE_GM_v0_34.r)
##            - code sample.NODEGM
##            - code predict.NODEGM
## 24-02-2022 - created v0_16
##            - cleaned code
##            - code storage of chain/ensemble
##            - clean code (esp. dependences, limit dependence on global variables)
## 03-03-2022 - created v0_17
##            - simplified preparation of data
## 16-03-2022 - created v0_18
##            - separate fitting script from visualisation script
##            - update fitting and plotting function with new repository
## 17-03-2022 - created v0_19
##            - impemented relative contributions
##            - implemented diamond plot
## 18-03-2022 - extra visualisations (diamond plots, phase space, ...)
##            - improve figures
## 20-03-2022 - created v0_20
##            - cleaned code by separating observation and process model
## 25-03-2022 - created v0_21
##            - implemented DEMC sampling of process model
##            - moved DEMC to annex file
## 28-03-2022 - created v0_22
##            - implemeted cross validation

## next:
##            - allow for different values of priors for each time series
##            - fix number of parameters in model
##            - allow to switch activation functions
##            - implement three-fold cross validation

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_50.r")

## load data
TS = read.table("data/TS_RPS.csv",sep=";",header=T)
s  = seq(0,600,10)
TS = TS[s,]
TS[,1] = TS[,1]/10
TS[,-1] = TS[,-1]/1000

# TS = read.table("data/TS_HL.csv",sep=";",header=T)
# TS[,1] = (TS[,1]-min(TS[,1]))
# TS[,-1] = TS[,-1]/10

## initiate time series 
t     = TS[,1]                                    # vector of time steps
Y     = TS[,-1]                                   # matrix containing time series

## output file
pdf("results.pdf")

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model 

## parameters
N_o    	 = 3      	       # number of ensemble elements
W_o    	 = 10    	       # number of nodes in observation model
sd1_o  	 = 0.1   	       # sd in likelihood of residuals of observation model
sd2_o  	 = 0.1    	       # sd in priors of observation model
alpha_i  = 1     	       # factor of interpolation

## interpolated time steps
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha_i)

## fit model_o 
# model_o = fit.NODEGM_o(t,t_,Y,W_o,sd1_o,sd2_o,N_o,log=T,logMar=T)
model_o = DFT.NODEGM_o(t,t_,Y,W_o,log=T)

## visualise
plot.NODEGM_o(t,t_,Y,model_o,col=rainbow(ncol(Y)))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model
 
## prepare data
Yhat_o     = NULL
ddt.Yhat_o = NULL
for(j in 1:(ncol(TS)-1))
{
    Yhat_o     = cbind(Yhat_o,    model_o[[j]]$MaP.Yhat_o)
    ddt.Yhat_o = cbind(ddt.Yhat_o,model_o[[j]]$MaP.ddt.Yhat_o)
}

## compute per-capita growth rate
rhat_o = 1/Yhat_o*ddt.Yhat_o
colnames(rhat_o) = colnames(Y)

## iterator 
for(i in 1:ncol(rhat_o))
{
	##
	message(paste("fitting: ",i,"/",ncol(rhat_o),sep=""))

	## variables
	X     = Yhat_o
	ddt.X = ddt.Yhat_o
	Y     = rhat_o[,i]
	
	## standardise predictive variables
	X_     = X
	mean_x = apply(X_,2,mean)
	sd_x   = apply(X_,2,sd)
	X_     = t((t(X_)-mean_x)/sd_x)
	
	## standardise response variable
	Y_     = Y
	mean_y = mean(Y_)
	sd_y   = sd(Y_)
	Y_     = (Y_-mean_y)/sd_y
			 
	## parameters
	N_p    	 = 3      	       # number of ensemble elements
	W_p    	 = 10     	       # number of nodes in process model
	sd1_p  	 = 0.1   	       # sd in likelihood of residuals of process model
	sd2_p  	 = 0.25  	       # sd in priors of process model
	
	## cross validation
	ensembleList = list()
	crossVal = NULL
	sd2_p_Vect = seq(0.1,0.5,0.05)
	for(k in 1:length(sd2_p_Vect))
	{
		##
		message(paste("crossval: ",k,"/",length(sd2_p_Vect),sep=""))

		## regularisation parameter
		sd2_p = sd2_p_Vect[k]
		
		## multi-fold cross validation
		nfolds    = 3
		folds     = rbind(c(1/4,2/4),c(3/4),c())
		crossVal_ = NULL
		for(u in 1:nfolds)
		{
			##
			message(paste("fold: ",u,"/",nfolds,sep=""))

			## split training/test
			s = sample(1:length(Y_),round(2/3*length(Y_))) # round(1/3*length(Y_)):round(3/3*length(Y_))
			#
			X_l = X_[s,]
			Y_l = Y_[s]
			#
			X_t = X_[-s,]
			Y_t = Y_[-s]

			## fit model_p
			ensemble_p = fit.model_p(X_l,Y_l,W_p,sd1_p,sd2_p,N_p)
			
			## model performance
			E.logLik_l  = mean(apply(ensemble_p,1,function(x)logLik_p(X_l,Y_l,x,sd1_p)))
			sd.logLik_l = sd(apply(ensemble_p,1,function(x)logLik_p(X_l,Y_l,x,sd1_p)))
			E.logLik_t  = mean(apply(ensemble_p,1,function(x)logLik_p(X_t,Y_t,x,sd1_p)))
			sd.logLik_t = sd(apply(ensemble_p,1,function(x)logLik_p(X_t,Y_t,x,sd1_p)))
		
			## cross validation matrix
			crossVal__ = c(sd2_p,E.logLik_l,sd.logLik_l,E.logLik_t,sd.logLik_t)
			crossVal_  = rbind(crossVal_,crossVal__)
		}

		## store ensemble
		ensembleList[[k]] = ensemble_p
		crossVal_ = apply(crossVal_,2,mean)
		crossVal = rbind(crossVal,crossVal_)
	}
	
	## visualise cross validation
	x = crossVal[,1]
	y = crossVal[,c(2,4)]
	plot(x,rep(0,length(x)),ylim=c(min(y)-0.1*abs(min(y)),max(y)+0.1*abs(max(y))),type="l",lty=3,xlab="sd",ylab="error")
	#
	y = crossVal[,2]; col="blue"
	sd.y = crossVal[,3]
	lines(x,y,col=col)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
	y = crossVal[,4]; col="red"
	sd.y = crossVal[,5]
	lines(x,y,col=col)
	polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
	#
	legend("topright",legend=c("training","testing"),col=c("blue","red"),lty=1,bty="n")
	
	## best regularisation
	k = which.max(crossVal[,4])
	ensemble_p = ensembleList[[k]]
	
	## predictions
	Yhat_p         = predict.model_p(X_,ensemble_p)
	ddx.Yhat_p     = ddx.predict.model_p(X_,ensemble_p)
	Geber_p        = lapply(ddx.Yhat_p,function(x) x * ddt.X)
	names(Geber_p) = c("E.Geber_p","q05.Geber_p","q95.Geber_p")
	
	## de-standardise predictions
	E.Yhat_p       = mean_y + sd_y * Yhat_p[["E.Yhat_p"]]
	q05.Yhat_p     = mean_y + sd_y * Yhat_p[["q05.Yhat_p"]]
	q95.Yhat_p     = mean_y + sd_y * Yhat_p[["q95.Yhat_p"]]
	#
	E.ddx.Yhat_p   = sd_y * t(1/sd_x * t(ddx.Yhat_p[["E.ddx.Yhat_p"]]))
	q05.ddx.Yhat_p = sd_y * t(1/sd_x * t(ddx.Yhat_p[["q05.ddx.Yhat_p"]]))
	q95.ddx.Yhat_p = sd_y * t(1/sd_x * t(ddx.Yhat_p[["q95.ddx.Yhat_p"]]))
	#
	E.Geber_p      = sd_y * t(1/sd_x * t(Geber_p[["E.Geber_p"]]))
	q05.Geber_p    = sd_y * t(1/sd_x * t(Geber_p[["q05.Geber_p"]]))
	q95.Geber_p    = sd_y * t(1/sd_x * t(Geber_p[["q95.Geber_p"]]))
	
	## figure
	xlab=c("","","Time");ylab=c("Dynamics","Effects","Contributions");index=NULL;legend=NULL;col="red"
	par(mfrow=c(3,1),mar=c(4,4,0,0),oma=c(1,1,1,1),cex.lab=1.25)
	#
	if (is.null(index)) index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
	N = ncol(E.ddx.Yhat_p)
	# 
	## dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(Y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=ylab[1])
	points(t_,Y) #,pch=16,col=adjustcolor("black",0.5))
	polygon(c(t_,rev(t_)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col,alpha=0.5),border=NA)
	lines(t_,E.Yhat_p,col=adjustcolor(col,alpha=0.75))
	if(!is.null(index)) legend("topright",legend=index[1],bty="n",cex=1)
	#
	## effects
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=ylab[2])
	for(j in 1:N) lines(t_,E.ddx.Yhat_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[2],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft",legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
	#
	## Geber 
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=ylab[3])
	for(j in 1:N) lines(t_, E.Geber_p[,j],col=rainbow(N,alpha=0.75)[j])
	for(j in 1:N) polygon(c(t_,rev(t_)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=rainbow(N,alpha=0.2)[j],border=NA)
	if(!is.null(index))  legend("topright",legend=index[3],bty="n",cex=1)
	if(!is.null(legend)) legend("topleft",legend=legend,lty=1,col=rainbow(N,alpha=0.75),bty="n")
	#
	par(mfrow=c(1,1),mar=c(4,4,3,3))
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
