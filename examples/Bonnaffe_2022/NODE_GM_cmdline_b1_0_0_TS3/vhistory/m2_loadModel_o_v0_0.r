#######################
## m2_loadModel_o.r ##
#######################

## goal: load observation model

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 06-04-2022 - created v0_0

## next:
##            - evaluate predictions for all Omega_iku
##            - allow modularity in process model
##            - fix number of parameters issue
##            - parallelise code

##############
## INITIATE ##
##############

## goal: initiate the NODE 

## load NODE functions
source("f_NODE_GM_v0_54.r")
source("m1_loadData_o.r")

## parameters
N = ncol(TS)-1

## parameters observation model
K_o      = 30                  # number of ensemble elements
W_o      = 10                  # number of nodes in observation model
N_o      = W_o * 3
sd1_o    = 0.1                 # sd in likelihood of residuals of observation model
sd2_o    = 0.1                 # sd in priors of observation model

#
###
