######################
## m6_loadModel_p.r ##
######################

## goal: load process model 

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0
##            - allow modularity in process model (create new f_p function in loadModel script)
##            - fix number of parameters issue
## 14-04-2022 - created v0_1
##            - fitting of process model is done on a different ensemble member of Yhat_o at every iteration
## 15-04-2022 - simplified code

##############
## INITIATE ##
##############

## goal: initiate the process model

## imports
source("f_NODE_GM_v0_56.r")
source("m5_loadData_p.r")

## parameters process model
K_p   = 10
W_p   = 10
N_p   = 2 * W_p * (2+N)
sd1_p = 0.1
sd2_p = c(0.3,0.3,0.3) # rep(0.1,N) 

#
###
