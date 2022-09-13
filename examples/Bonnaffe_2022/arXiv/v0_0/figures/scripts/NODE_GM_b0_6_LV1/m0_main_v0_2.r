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

## goal: train process model 
source("m7_crossVal_p.r")
source("m8_plotCrossVal_p.r")

## goal: visualise process model results
source("m9_trainModel_p.r")
source("m10_plotModel_p.r")

#
###
