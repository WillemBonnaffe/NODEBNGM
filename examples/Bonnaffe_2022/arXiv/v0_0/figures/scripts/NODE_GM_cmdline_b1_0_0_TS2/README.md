############
## README ##
############

## update log:
## 28-03-2022 - create b0_6
##            - implemented cross validation
## 12-04-2022 - created b0_9
##            - this is the best version
##            - improve replicability of results though
## 15-04-2022 - created b0_12
##            - implemented dataloaders
##            - implemented dataloaders with random selection of observation data
##            - implemented bootstrapping of observation data in the cross validation and training process
## 25-04-2022 - created b0_13
##            - simplify and clean code
##              simplified code drastically (1000 -> 500 lines in function file)
##              the code of the observation model is quite clean
##              the code of the process model is still a bit messy (esp crossval)
##            - cleaned module names
##              harmonise the file names in the repository (quite messy)
##            - fixate sd and mean of data in dataloader (think about what it means to standardise each response and predictive variable coming for the observation distribution)
##            - fix error in train model (sd and mean are not matching the predictions to the right data)
##              if too much sampling is required for crossvalidation with random dataloader then can use expectation of data 
##              for this make sure that many samples are taken for the observation model (to ensure repeatability of mean)
##            - improve management of output folder (rn it is called by the dataloader which is not good)

## to do next:
##            - think about method
##            - test further the repeatability of crossvalidation and results
##            - can use also just b0_9 and increase the size of the samples taken for interpolation to increase repeatability of results
##            - investigate nonlinearity in LV4

## method:
## 1. Observation model - works well but generate more samples
## 2. Process model 
##    * crossvalidation on expectation of observation model (otherwise requires too many 
##      samples)
##    * train process model with uncertainty of the observation model

## goals:
## x strip code to essential elements (remove any extra functions)
## x fix errors in code
## * re-think the structure of the repository to allow modularity and parallelisation
## * ensure repeatability of model by re-running approach a second time (esp for cross
##   validation on expectation of observation model)
## * assess capacity of model to recover linear response for linear model
## * assess capacity of model to recover nonlinear response for nonlinear model

#
###
