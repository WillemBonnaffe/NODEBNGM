#####################
## m1_loadData_o.r ##
#####################

## goal: load time series and prepare for NODE analysis 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 06-04-2022 - created v0_0
## 15-04-2022 - renamed file m1_loadData_o.r
##            - created v0_1
##            - implemented peparation of data in the dataloader

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## TS 3DLV ##
## load data
TS = read.table("data/TS_3DLV.csv",sep=";",header=T)
TS = TS[20:50,]
# TS = TS[1:30,]
#
## make out directory
pathToOut = "out/"
system(paste("mkdir",pathToOut))

## parameters
alpha_i = 1

#
###

##################
## PREPARE DATA ##
##################

## goal: prepare data for training of observation model

## data specs
N       = ncol(TS) - 1
n       = nrow(TS)

## prepare data
X_o = TS[,1]
Y_o = TS[,-1]

## predictive and response variable
t   = X_o 
Y   = Y_o
nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)

## standardise time steps
t_0 = min(t)
t_f = max(t)
dt  = diff(t[1:2])
t_  = 1:length(t)
nt_ = seq(min(t_),max(t_),(t_[2]-t_[1])/alpha_i)

## standardise data
Y_     = log(Y)
mean_y = apply(Y_,2,mean)
sd_y   = apply(Y_,2,sd)
Y_     = t((t(Y_)-mean_y)/sd_y)

#
###
