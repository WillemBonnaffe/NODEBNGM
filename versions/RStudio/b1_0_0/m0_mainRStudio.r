######################
## m0_mainRStudio.r ##
######################

## goal: run all modules

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

#
###

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load data
TS = read.table("data/TS_1.csv",sep=";",header=T)
for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.005)] = 0.005}

## make out directory
pathToOut = "out/"
system(paste("mkdir",pathToOut))

## load NODE functions
source("f_NODE_GM_Rstudio.r")

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## parameters
K_o     = 100
W_o     = c(100,100,100)
N_o     = W_o*3
rho     = 1
alpha_i = 1

## train 
results_o = trainModel_o(TS,alpha_i,N_o,K_o,rho)
attach(results_o,warn.conflicts=F)

## visualise
pdf(paste(pathToOut,"/fig_predictions_o.pdf",sep=""))
plotModel_o(TS,alpha_i,Yhat_o,ddt.Yhat_o)
dev.off()

## save results
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## parameters
K_p   = 10
W_p   = rep(10,N)
N_p   = 2 * W_p * (2+N)
sd1_p = 0.1
sd2_p = list(c(rep(1.0,N_p[1]/2),rep(.15,N_p[1]/2)),
             c(rep(1.0,N_p[2]/2),rep(.01,N_p[2]/2)),
             c(rep(1.0,N_p[3]/2),rep(.075,N_p[3]/2)))

## train 
results_p = trainModel_p(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p)
attach(results_p,warn.conflicts=F)

## visualise 
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
dev.off()

## store results  
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))

#
###

###################################
## CROSSVALIDATION PROCESS MODEL ##
###################################

## parameters
K_p             = 3
folds           = list(c(1/2,1)) # list(c(0,1/3),c(1/3,2/3),c(2/3,3/3))
crossValParVect = seq(0.01,0.5,0.01)

## cross validation
resultsCrossVal_p = crossVal_p(TS,alpha_i,Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,folds,crossValParVect)

## visualise
pdf(paste(pathToOut,"/fig_crossVal_p.pdf",sep=""))
plotCrossVal_p(resultsCrossVal_p)
dev.off()

## store results
save(resultsCrossVal_p,file=paste(pathToOut,"/","crossVal_p.RData"   ,sep=""))

#
###
