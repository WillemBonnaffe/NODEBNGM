###############
## m0_main.r ##
###############

## goal: run all modules

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 09-06-2022 - created v0_0

## goal: train observation model
source("m3_trainModel_o.r")
source("m4_plotModel_o.r")

## goal: cross validation on process model
source("m9_crossVal_p.r")
source("m10_plotCrossVal_p.r")

## goal: train process model 
source("m7_trainModel_p.r")
source("m8_plotModel_p.r")

#
###
