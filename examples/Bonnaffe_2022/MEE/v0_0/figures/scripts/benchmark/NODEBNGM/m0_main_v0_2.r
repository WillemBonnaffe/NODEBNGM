###############
## m0_main.r ##
###############

## goal: run all modules

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0
## 15-04-2022 - created v0_1
##            - updated dependencies
##            - fix number of parameters issue
##            - allow modularity in process model => can specify custom f_p in loadModel script
## 25-04-2022 - created v0_2
##            - updated name of files
##            - included cross validation

## next:
##            - parallelise code

## goal: train observation model
source("m3_trainModel_o.r")

## goal: visualise observation model results
source("m4_plotModel_o.r")

# ## goal: train process model 
source("m7_trainModel_p.r")
# source("m7b_predictModel_p.r")

## goal: visualise process model results
source("m8_plotModel_p.r")

## goal: cross validation
source("m9_crossVal_p.r")
source("m10_plotCrossVal_p.r")

#
###
