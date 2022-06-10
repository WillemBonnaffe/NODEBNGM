###############
## m0_main.r ##
###############

## goal: run all modules

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 10-04-2022 - created v0_0

## module m1
## goal: load time series data
source("m1_loadData.r")

## module m2
## goal: load observation model 
source("m2_loadModel_o.r")

## module m3
## goal: train observation model
source("m3_trainModel_o.r")

## module m4
## goal: visualise observation model results
source("m4_visModel_o.r")

## module m5
## goal: load process model 
source("m5_loadModel_p.r")

## module m6
## goal: train process model 
source("m6_trainModel_p.r")

## module m7
## goal: visualise process model results
source("m7_visModel_p.r")

#
###
