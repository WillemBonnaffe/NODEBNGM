#####
## ## 
#####


data(block_3sp)
lib <- c(1, 100)
pred <- c(101, 198)

cols <- c("x_t", "y_t", "z_t")
target <- "x_t"

block_smap_output <- block_lnlp(block_3sp, lib = lib, pred = pred, columns = cols, 
                                target_column = target, method = "s-map", theta = 2, stats_only = FALSE, 
                                first_column_time = TRUE, save_smap_coefficients = TRUE, silent = TRUE)

smap_coeffs <- block_smap_output$smap_coefficients[[1]]
str(smap_coeffs)

predictions <- block_smap_output$model_output[[1]]
t <- predictions$time

plot(t, predictions$obs, type = "l", col = "black", ylab = "x", xlab = "")
lines(t, predictions$pred, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")
#
par(mfrow=c(3,1))
plot(t, smap_coeffs[, 2], type = "l", col = "red", ylab = "effect of x", xlab = "")
plot(t, smap_coeffs[, 3], type = "l", col = "blue", ylab = "effect of y", xlab = "")
plot(t, smap_coeffs[, 4], type = "l", col = "magenta", ylab = "effect of z", xlab = "")
par(mfrow=c(1,1))

##############
## INITIATE ##
##############

## 
library(rEDM)

## load data
data(block_3sp)
lib <- c(1, 100)
pred <- c(101, 120)
cols <- c("x_t", "y_t", "z_t")
target <- "x_t"

## load data
TS = as.matrix(read.table("data/TS_3DLV.csv",sep=";",header=T))
TS = TS# [20:50,]
# TS[,-1] = log(TS[,-1])# [20:50,]
TS[,1] = TS[,1]-min(TS[,1])

## format for smaps
TS = data.frame(TS)
lib <- c(1, 101)
pred <- c(1, 101)
cols <- c("G", "B", "R")
target <- "B"

## load ground truth
load("data/GT_3DLV.RData")
ddx.rhat_true_t = LV_GT$t_true
ddx.rhat_true = LV_GT$ddx.rhat_true
i = 2

#
###

##########
## MAIN ##
##########

## S-map
block_smap_output <- block_lnlp(TS, 
                                lib = lib, 
                                pred = pred, 
                                columns = cols, 
                                target_column = target, 
                                method = "s-map", 
                                theta = 2, 
                                stats_only = FALSE, 
                                first_column_time = TRUE, 
                                save_smap_coefficients = TRUE, 
                                silent = TRUE)
smap_coeffs <- block_smap_output$smap_coefficients[[1]]
str(smap_coeffs)

## prediction
predictions <- block_smap_output$model_output[[1]]
t <- predictions$time

## figure
#
## predictions
plot(predictions$t, predictions$Observations, type = "l", col = "black", ylab = "x", xlab = "")
lines(predictions$t, predictions$Predictions, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")
#
## Jacobian
plot(1:10,xlim=c(0,max(smap_coeffs$t)),ylim=c(-1,1)*3,cex=0)
#
correction = 1/LV_GT$Yhat_true[,i]*LV_GT$ddt.Yhat_true[,i]*0
#
lines(smap_coeffs$t, smap_coeffs[, 3],   col = "green")
lines(ddx.rhat_true_t,LV_GT$ddx.rhat_true[,1+3*(i-1)], col = "green", lty=2)
# lines(ddx.rhat_true_t,LV_GT$Yhat_true[,i]*LV_GT$ddx.rhat_true[,1+3*(i-1)] + correction*(i==1), col = "green", lty=2)
lines(smap_coeffs$t, smap_coeffs[, 4],   col = "blue")
lines(ddx.rhat_true_t,LV_GT$ddx.rhat_true[,2+3*(i-1)], col = "blue", lty=2)
# lines(ddx.rhat_true_t,LV_GT$Yhat_true[,i]*LV_GT$ddx.rhat_true[,2+3*(i-1)] + correction*(i==2), col = "blue", lty=2)
lines(smap_coeffs$t, smap_coeffs[, 5],   col = "red")
lines(ddx.rhat_true_t,LV_GT$ddx.rhat_true[,3+3*(i-1)], col = "red", lty=2)
# lines(ddx.rhat_true_t,LV_GT$Yhat_true[,i]*LV_GT$ddx.rhat_true[,3+3*(i-1)] + correction*(i==3), col = "red", lty=2)

## dY/dt = r Y
## ddy(dY/dt) = ddy(r Y)
## ddy(dY/dt) = ddy(r) Y + r ddy(Y)
## ddy(r) = 1/Y (ddy(dY/dt) - r)
## ddy(r) = 1/Y (J - 1/Y dY/dt)

## substitution
## X = log(Y) => dX/dY = 1/Y
##
## derivative wtr t
## dX/dt = d(log(Y))/dt
## dX/dt = dY/dt 1/Y
## dY/dt = Y dX/dt
## dY/dt = exp(X) dX/dt
##
## computing the Jacobian
## d(dY/dt)/dY = d(dY/dt)/dX dX/dY (chain rule)
## d(dY/dt)/dY = d(exp(X) dX/dt)/dX 1/Y
## d(dY/dt)/dY = (exp(X) dX/dt + exp(X) d(dX/dt)/dX) 1/Y
## d(dY/dt)/dY = (exp(X)/Y dX/dt + exp(X)/Y d(dX/dt)/dX)
## d(dY/dt)/dY = dX/dt + d(dX/dt)/dX
## J_Y = dX/dt + J_X
## J_X = J_Y - r

#
###