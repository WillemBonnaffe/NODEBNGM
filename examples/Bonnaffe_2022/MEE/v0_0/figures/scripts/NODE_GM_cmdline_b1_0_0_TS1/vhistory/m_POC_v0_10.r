###########
## m_POC ##
###########

## goal: perform assessment of quality of inference of NODE by gradient matching 

## update:
## 03-03-2022 - created v0_0 from m_NODE_GM_v0_17.r
## 07-03-2022 - created v0_1
##            - re-introduced scaling of variables to help with the fitting 
##            - created v0_2
##            - tested DEMCpp and DEMCO
## 08-03-2022 - created v0_3
##            - coded anchor sampling
##            - added specific priors
##            - re-introduced FFT as benchmark
## 10-03-2022 - created v0_5
##            - cleaned code (removed unused sections)
## 11-03-2022 - created v0_6
##            - implemented FFT and DFT
##            - added tri-trophic time series 
## 14-03-2022 - created v0_7
##            - individual scaling descaling
##            - compressed and cleaned code
## 14-03-2022 - created v0_8
##            - tweaked approach to work with RPS time series (simply changed units in time series)
##            - move functions to f_ file
## 15-03-2022 - created v0_9
##            - implemented main functions created from moving parts of the code into functions
##            - created v0_10
##            - compressed code
##            - modularise code (use compact functions, handles the scaling)
##            - separate function logPost and logMar

## next:
##            - clean and compress code 
##            - create plotting functions
##            - create markdown
##            - create compact code

##############
## INITIATE ##
##############

### load NODE functions
source("f_NODE_GM_v0_41.r")

## 
library(deSolve)

## TS to consider (1,2,3)
TS_id = 1

## output file
pdf(paste("results_TS_",TS_id,".pdf",sep=""))

#
###

############################
## ARTIFICIAL TIME SERIES ##
############################

if (TS_id == 1)
{
	#############
	## MODEL 1 ##
	#############
	
	## goal: 
	
	## run ode
	t_true     = seq(0,20,0.1)
	Y_true_0   = c(1.0,1.0)
	p          = c(1.0,0.5,1.0)
	ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]*x[2],p[2]*x[1]*x[2]-p[3]*x[2]))
	ddx.r_true = function(t,x,p) list(c(0,-p[2],p[2],0))
	Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)
	
	## create TS 
	s = seq(0,20*10,5)
	TS = Yhat_true[s,]
	colnames(TS) = c("t","hare","lynx")
	write.table(TS,file="data/TS_LV_1.csv",sep=";",row.names=FALSE)
	
	## ground truth
	Yhat_true     = Yhat_true[,2:3]
	ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
	ddx.r_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))
	
	## visualise
	plot(TS[,1],TS[,2],col="red",type="l")
	lines(TS[,1],TS[,3],col="blue")

	#
	###
}
 
if (TS_id == 2)
{
	#############
	## MODEL 2 ##
	#############
	
	## goal: simple LV with no interaction
	
	## run ode
	t_true     = seq(0,20,0.1)
	Y_true_0   = c(1.0,1.0)
	p          = c(1.0,1.5,0.25,1.0)
	ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]/(1+p[2]*p[3]*x[1])*x[2],p[2]*x[1]/(1+p[2]*p[3]*x[1])*x[2]-p[4]*x[2]))
	ddx.r_true = function(t,x,p) list(c(+p[2]*x[2]*p[2]*p[3]/(1+p[2]*p[3]*x[1])^2,
									    -p[2]/(1+p[2]*p[3]*x[1]),
	                                     p[2]/(1+p[2]*p[3]*x[1])^2,
	                                     0))
	Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)
	
	## create TS
	s = seq(1,20*10,5)
	TS = Yhat_true[s,]
	colnames(TS) = c("t","hare","lynx")
	write.table(TS,file="data/TS_LV_2.csv",sep=";",row.names=FALSE)
	
	## ground truth
	Yhat_true     = Yhat_true[,2:3]
	ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
	ddx.r_true = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))
	
	## visualise
	plot(TS[,1],TS[,2],col="red",type="l")
	lines(TS[,1],TS[,3],col="blue")
	
	#
	###
}

if (TS_id == 3)
{
	#############
	## MODEL 3 ##
	#############
	
	## goal: simple LV with no interaction
	
	## run ode
	t_true     = seq(0,20,0.1)
	Y_true_0   = c(1.0,1.0,0.1)
	p          = c(1.0,0.5,0.25,0.4,1.0,1.0)
	ddt.Y_true = function(t,x,p) list(c(p[1]*x[1]-p[2]*x[1]*x[2]-p[3]*x[1]*x[3],
									              p[2]*x[1]*x[2]-p[4]*x[2]*x[3]-p[5]*x[2],
										          p[3]*x[1]*x[3]+p[4]*x[2]*x[3]-p[6]*x[3]))
	ddx.r_true = function(t,x,p) list(c(0,-p[2],-p[3],p[2],0,-p[4],p[3],p[4],0))
	Yhat_true  = ode(y=Y_true_0,times=t_true,func=ddt.Y_true,parms=p)
	
	## create TS 
	s = seq(0,20*10,5)
	TS = Yhat_true[s,]
	colnames(TS) = c("t","R","G","B")
	write.table(TS,file="data/TS_LV_3.csv",sep=";",row.names=FALSE)
	
	## ground truth
	Yhat_true     = Yhat_true[,2:4]
	ddt.Yhat_true = t(apply(Yhat_true,1,function(x)unlist(ddt.Y_true(NULL,x,p))))
	ddx.r_true    = t(apply(Yhat_true,1,function(x)unlist(ddx.r_true(NULL,x,p))))
	
	## visualise
	plot(TS[,1],TS[,2],col="red",type="l",ylim=c(0,max(TS[,-1])))
	lines(TS[,1],TS[,3],col="green")
	lines(TS[,1],TS[,4],col="blue")
	
	#
	###
}

#
###

##########################
## INITIATE TIME SERIES ##
##########################

## goal: initiate the time series

## load data
# TS = read.table("data/TS_LV_3.csv",sep=";",header=T)

#
###
 
###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model with DFT and FFT

## parameters 
K     = 10
alpha = 2	

## initiate time series
t     = TS[,1]
t_    = seq(min(t),max(t),(t[2]-t[1])/alpha)

## for each column in TS
pred_o_DFT = list()
for(i in 1:(ncol(TS)-1))
{
	## response variable
	Y = TS[,i+1]

	## DFT
	DFT_            = fit.DFT(Y,t,t_,K,log=T)
	Yhat_o_DFT      = DFT_[,1]
	ddt.Yhat_o_DFT  = DFT_[,2]
	pred_o_DFT[[i]] = list("Yhat_o_DFT"      = Yhat_o_DFT,
							"ddt.Yhat_o_DFT" = ddt.Yhat_o_DFT)
	
	## visualise
	par(mfrow=c(2,1))
	#
	plot(t_,rep(0,length(t_)),ylim=c(0,1)*max(abs(TS[,i+1]))*1.5,type="l",lty=3)
	points(t,TS[,i+1])
	lines(t_,Yhat_o_DFT,col="blue",lty=1)
	#
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(ddt.Yhat_o_DFT))*1.5,type="l",lty=3)
	lines(t_,ddt.Yhat_o_DFT,col="blue",lty=1)
	lines(t_true,ddt.Yhat_true[,i],type="l",lty=2,ylim=c(-1,1)*3)
	#
	par(mfrow=c(1,1))
}

#
###

###########################
## FIT OBSERVATION MODEL ## 
###########################

## goal: fit observation model

## parameters
N_e    = 3    					                    # number of ensemble elements
W_o    = 10    					                    # number of nodes in observation model
sd1_o  = 0.01   					                # sd in likelihood of residuals of observation model
sd2_o  = c(rep(0.2,W_o),rep(0.1,W_o),rep(0.1,W_o))	# sd in priors of observation model
alpha  = 1      	   								# factor of interpolation

## initiate time steps
t  = TS[,1] 		 	                            # vector of time steps
t_ = seq(min(t),max(t),(t[2]-t[1])/alpha)  	        # vector of interpolated time steps

## for each column in TS
pred_o = list()
for (i in 1:(ncol(TS)-1))
{
	## 
	message(paste("fitting: ",i,"/",ncol(TS)-1,sep=""))

	## response variable
	Y = TS[,i+1]
	
	## fit observation model
	Omega_o_chain  = fit.model_o(t,Y,W_o,sd1_o,sd2_o,N_e,log=T,logMar=T)
	pred_o[[i]]    = predict.model_o(t_,Omega_o_chain,log=T)
	attach(pred_o[[i]],warn.conflicts=F)

	## visualise density and dynamics
	par(mfrow=c(2,1))
	#
	## plot density
	plot(t_,rep(0,length(t_)),ylim=c(0,1)*max(Y)*1.5,type="l",lty=3)
	points(t,Y)
	polygon(x=c(t_,rev(t_)),y=c(q05.Yhat_o,rev(q95.Yhat_o)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.Yhat_o,col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_,MaP.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i],lty=2)
	lines(t_true,Yhat_true[,i],lty=2)
	#
	## plot dynamics
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(MaP.ddt.Yhat_o))*1.5,type="l",lty=3)
	polygon(x=c(t_,rev(t_)),y=c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.ddt.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_,MaP.ddt.Yhat_o,type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[i],lty=2)
	lines(t_true,ddt.Yhat_true[,i],lty=2)
	#
	par(mfrow=c(1,1))
	
	## visualise periodic components
	Omega_o_MaP = Omega_o_chain[which.max(Omega_o_chain[,1]),-c(1:3)]
	Yhat_o_     = t(apply(matrix(Omega_o_MaP,ncol=3),1,function(x) f_o.eval(t_,x)))
	sum.Yhat_o_ = apply(Yhat_o_,2,sum)
	plot(t_,sum.Yhat_o_,ylim=c(min(Yhat_o_),max(Yhat_o_)),cex=0)
	for(j in 1:W_o)
	{
		lines(t_,Yhat_o_[j,],col=rainbow(W_o)[j])
	}
	 
	## visualise prior and posterior distributions
	par(mfrow=c(3,1))
	for(j in 1:3)
	{
		chain_ = as.vector(Omega_o_chain[,-1][,1:W_o+W_o*(j-1)])
		sd2_o_ = sd2_o[1+W_o*(j-1)]
		x      = density(chain_)$x
		y      = density(chain_)$y
		pd     = 1/sqrt(2*pi*sd2_o_^2)*exp(-(x^2)/(2*sd2_o_^2))
		plot(x,pd,type="l",lty=2,ylim=c(0,max(c(y,pd))))
		lines(x,y)
	}
	par(mfrow=c(1,1))

}

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## parameters
N_e    	 = 3     					# number of ensemble elements
W_p    	 = 10     					# number of nodes in process model
sd1_p  	 = 0.01   	   				# sd in likelihood of residuals of process model
sd2_p  	 = 0.05  	   				# sd in priors of process model

## prepare data
Yhat_o     = NULL 
ddt.Yhat_o = NULL
for(i in 1:(ncol(TS)-1))
{
	Yhat_o     = cbind(Yhat_o,    pred_o[[i]]$MaP.Yhat_o)
	ddt.Yhat_o = cbind(ddt.Yhat_o,pred_o[[i]]$MaP.ddt.Yhat_o)
}
#
## or use DFT
# Yhat_o     = Yhat_o_DFT 
# ddt.Yhat_o = ddt.Yhat_o_DFT 
#
## or ground truth
# Yhat_o     = ode(y=Y_true_0,times=t_,func=ddt.Y_true,parms=p)[,2:3]
# ddt.Yhat_o = t(apply(Yhat_o,1,function(x)unlist(ddt.Y_true(NULL,x,p))))

## compute per-capita growth rate
r_o = 1/Yhat_o*ddt.Yhat_o

pred_p = list()
Geber_p = list()
for(i in 1:(ncol(TS)-1))
{
	## 
	message(paste("fitting: ",i,"/",ncol(TS)-1,sep=""))

	## data
	X      = Yhat_o
	ddt.X  = ddt.Yhat_o
	Y      = r_o[,i]
	
	## predictions
	Omega_p_chain = fit.model_p(X,Y,W_p,sd1_p,sd2_p,N_e)
	pred_p[[i]]   = predict.model_p(X,Omega_p_chain)
	Geber_p[[i]]  = Geber.model_p(X,ddt.X,Omega_p_chain)
	attach(pred_p[[i]],warn.conflicts=F)
	attach(Geber_p[[i]],warn.conflicts=F)

	## visualise predictions
	par(mfrow=c(3,1))
	#
	## plot y
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(r_o))*1.5,type="l",lty=3)
	points(t_,Y)
	polygon(x=c(t_,rev(t_)),y=c(q05.Yhat_p,rev(q95.Yhat_p)),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[i])
	lines(t_,E.Yhat_p,col=rainbow((ncol(TS)-1),alpha=0.5)[i])
	lines(t_true,1/Yhat_true[,i]*ddt.Yhat_true[,i],lty=2)
	#
	## plot derivative
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3)
	for(j in 1:(ncol(TS)-1))
	{
		polygon(x=c(t_,rev(t_)),y=c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[j])
		lines(t_,E.ddx.Yhat_p[,j],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j])
		lines(t_true,ddx.r_true[,j+(i-1)*(ncol(TS)-1)],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j],lty=2)
	}
	#
	## plot Geber
	plot(t_,rep(0,length(t_)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3)
	for(j in 1:(ncol(TS)-1))
	{
		polygon(x=c(t_,rev(t_)),y=c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),border=NA,col=rainbow((ncol(TS)-1),alpha=0.25)[j])
		lines(t_,E.Geber_p[,j],type="l",col=rainbow((ncol(TS)-1),alpha=0.5)[j])
	}
	#
	par(mfrow=c(1,1))
	
	## visualise prior and posterior distributions
	par(mfrow=c(3,1))
	for(j in 1:3)
	{
		tmp = as.vector(Omega_p_chain[,-1][,1:W_p+W_p*(j-1)])
		# sd2_p_ = sd2_[1+W_p[i]*(j-1)]
		x = density(tmp)$x
		y = density(tmp)$y
		pd = 1/sqrt(2*pi*sd2_p^2)*exp(-(x^2)/(2*sd2_p^2))
		plot(x,pd,type="l",lty=2,ylim=c(0,max(c(y,pd))))
		lines(x,y)
	}
	par(mfrow=c(1,1))
	
	# ## visualise NODE approx
	# par(mfrow=c(2,2))
	# #
	# func = function(x)
	# {
	# 	x = mean(X)  + sd(X) * x 
	# 	r = 1/x[i]*unlist(ddt.Y_true(NULL,x,p))[i]
	# 	return(r)
	# }
	# .plotMat(rbind(c(-3,3),c(-3,3)),func,maxAbsMinMax=5)
	# #
	# func = function(x) predict.model_p(t(x),Omega_p_chain)$E.Yhat_p
	# .plotMat(rbind(c(-3,3),c(-3,3)),func,maxAbsMinMax=5)
	# #
	# par(mfrow=c(1,1))
}

#
###

dev.off()
