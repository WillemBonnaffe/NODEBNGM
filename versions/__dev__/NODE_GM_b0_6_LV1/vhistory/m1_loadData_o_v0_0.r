###################
## m1_loadData_o.r ##
###################

## goal: load time series and prepare for NODE analysis 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 06-04-2022 - created v0_0
## 15-04-2022 - renamed file m1_loadData_o.r

## next:
##            - evaluate predictions for all Omega_iku
##            - fix number of parameters issue
##            - allow modularity in process model
##            - parallelise code

##############
## INITIATE ##
##############

## goal: initiate the NODE 

# ## algae flagellate rotifer time series ##
# ## load data
# TS = read.table("data/TS_1.csv",sep=";",header=T)
# #
# ## avoid nan errors when logging
# for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.005)] = 0.005}
# #
# ## make out directory
# pathToOut = "out/test/"
# system(paste("mkdir",pathToOut))

## rock paper scissor time series ##
## load data
TS = read.table("data/TS_RPS.csv",sep=";",header=T)
s  = seq(0,600,10)
TS = TS[s,]
TS[,1] = TS[,1]/10
TS[,-1] = TS[,-1]/1000
#
## make out directory
pathToOut = "out/RPS/"
system(paste("mkdir",pathToOut))

# ## short hare-lynx time series ##
# ## load data
# TS = read.table("data/TS_HL.csv",sep=";",header=T)
# s  = 1:60
# TS = TS[-s,]
# TS[,1] = (TS[,1]-min(TS[,1]))
# TS[,-1] = TS[,-1]/20
# #
# ## make out directory
# pathToOut = "out/HL_short/"
# system(paste("mkdir",pathToOut))

# ## long hare-lynx time series ##
# ## load data
# TS = read.table("data/TS_HL.csv",sep=";",header=T)
# s  = 1:30
# TS = TS[-s,]
# TS[,1] = (TS[,1]-min(TS[,1]))
# TS[,-1] = TS[,-1]/20
# #
# ## make out directory
# pathToOut = "out/HL_long/"
# system(paste("mkdir",pathToOut))

#
###
