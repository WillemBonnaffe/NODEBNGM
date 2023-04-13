#######################
## m0_main_RStudio.r ##
#######################

## goal: 
## - use Bayesian neural gradient matching (BNGM) to fit neural ordinary differential equation model (NODE) to time series 
## - analyse NODE to derive interactions and contributions between variables (sensitivities, Geber method)

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## update log:
## 01-12-2022 - created v0_0
## 12-04-2023 - implemented more than two fold cross validation
## 13-04-2023 - simplified code for use as template

#
###

##############
## INITIATE ##
##############

## goal: load data, functions

## load NODEBNGM functions
source("f_NODEBNGM_Rstudio.r")

## load data
TS = read.table("data/TS_Ushio.csv",sep=",",header=T)

## extract time steps and columns of interest
selected_time_steps = 50:150
selected_columns  = c(
  "time_step",
  # "surf.t",
  "bot.t",
  # "Aurelia.sp",                  # removed for example
  # "Engraulis.japonicus",         # too sparse
  # "Plotosus.lineatus",           # too sparse
  # "Sebastes.inermis",            # removed for example
  "Trachurus.japonicus",
  # "Girella.punctata",            # removed for example
  "Pseudolabrus.sieboldi"#,
  # "Halichoeres.poecilopterus",   # removed for example
  # "Halichoeres.tenuispinnis",    # removed for example
  # "Chaenogobius.gulosus",        # too sparse
  # "Pterogobius.zonoleucus",      # removed for example
  # "Tridentiger.trigonocephalus", # removed for example
  # "Siganus.fuscescens",          # too sparse
  # "Sphyraena.pinguis",           # too sparse
  # "Rudarius.ercodes"             # removed for example
)
TS = TS[selected_time_steps,]
TS = TS[,selected_columns]

## shorten column names
column_names =  c("time_step",
                    "bot.t",
                    # "Aurel.sp",
                    # "S.inerm.",
                    "T.japon.",
                    # "G.punct.",
                    "P.siebo."#,
                    # "H.poeci.",
                    # "H.tenui.",
                    # "P.zonol.",
                    # "T.trigo.",
                    # "S.pingu.",
                    # "R.ercod."
                  )
colnames(TS) = column_names

## normalise time series
TS[,-1] = apply(TS[,-1],2,function(x)(x-min(x))/(max(x)-min(x))*10)

## set 0s to small value to avoid NAs
for(i in 2:ncol(TS)){TS[,i][which(TS[,i]<0.005)] = 0.005}

## make output directory
pathToOut = "out"
system(paste("mkdir",pathToOut))

## visualise time series
pdf(paste(pathToOut,"/fig_time_series.pdf",sep=""))
par(mfrow=c(3,4))
for(i in 2:ncol(TS))
{
  plot(TS[,1],TS[,i],type="l",xlab="Time step",ylab="Relative density",bty="n",main=colnames(TS)[i])
}
par(mfrow=c(1,1))
dev.off()

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## goal: fit observation model (i.e. interpolate each variable in time series and compute interpolated temporal derivative to approximate temporal dynamics)

## parameters of observation model
N       = ncol(TS) - 1
K_o     = 3                  # number of ensemble elements
W_o     = rep(30,N)          # number of neurons in observation model, by default a single layer perceptron (equivalent to number of elements in Fourier series)
N_o     = W_o*3              # total number of parameters in observation model
rho     = 1                  # proportion of best samples to reject (to refine quality of fit if necessary)
alpha_i = 1                  # upsampling interpolation factor (2 double the number of points in interpolated time series)

## train observation model
runtime_o = system.time({
    model_o     = trainModel_o(TS,alpha_i,N_o,K_o,rho)
})[3]
Yhat_o      = model_o$Yhat_o
ddt.Yhat_o  = model_o$ddt.Yhat_o
Omega_o     = model_o$Omega_o

## visualise observation model fit
pdf(paste(pathToOut,"/fig_predictions_o.pdf",sep=""))
plotModel_o(TS,alpha_i,Yhat_o,ddt.Yhat_o)
dev.off()

## save results
save(Yhat_o,    file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
save(ddt.Yhat_o,file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
save(Omega_o,   file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))
write(runtime_o,   file=paste(pathToOut,"/","runtime_o.txt",sep=""))

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## goal: fit process model (i.e. explain the per-capita growth rate of the populations calculated as 1/Y*dY/dt as a function of the states Y(t))

## load model_o
load(file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
load(file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
load(file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

## parameters of process model
K_p   = 3                      # number of models to fit
W_p   = rep(10,N)              # number of neurons in single layer perceptron (SLP)
N_p   = 2 * W_p * (2+N)        # number of parameters in process model
sd1_p = 0.1                    # standard deviation of model likelihood
sd2_p = as.list(rep(0.1,N))   # standard deviation of prior distributions (second half concerns nonlinear functions)
         
## train process model
runtime_p = system.time({
    model_p    = trainModel_p(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,trainSplit=2/3)
})
Yhat_p     = model_p$Yhat_p     
ddx.Yhat_p = model_p$ddx.Yhat_p 
Geber_p    = model_p$Geber_p   
Omega_p    = model_p$Omega_p   

## remove drivers of temperature as not driven by variables
Yhat_p[[1]]     = Yhat_p[[1]]*0
ddx.Yhat_p[[1]] = ddx.Yhat_p[[1]]*0
Geber_p[[1]]    = Geber_p[[1]]*0
Omega_p[[1]]    = Omega_p[[1]]*0

## visualise process model
pdf(paste(pathToOut,"/fig_predictions_p.pdf",sep=""))
plotModel_p(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
dev.off()

## store results 
save(Yhat_p       ,file=paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
save(ddx.Yhat_p   ,file=paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
save(Geber_p      ,file=paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
save(Omega_p      ,file=paste(pathToOut,"/","Omega_p.RData"   ,sep=""))
write(runtime_p,   file=paste(pathToOut,"/","runtime_p.txt",   sep=""))

#
###

##############
## ANALYSIS ##
##############

## load results process model
load(paste(pathToOut,"/","Yhat_p.RData"    ,sep=""))
load(paste(pathToOut,"/","ddx.Yhat_p.RData",sep=""))
load(paste(pathToOut,"/","Geber_p.RData"   ,sep=""))
load(paste(pathToOut,"/","Omega_p.RData"   ,sep=""))

## compute Jacobian and contribution matrix
MSq = function(x) mean(x^2)
prop = function(x) x/sum(x)
J = t(matrix(unlist(lapply(ddx.Yhat_p,function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,mean))),ncol=ncol(TS)-1)) ## average across samples then average across time steps
C = t(matrix(unlist(lapply(Geber_p,   function(x)apply(matrix(apply(x,2,mean),nrow=nrow(TS),byrow=T),2,MSq))),ncol=ncol(TS)-1)) ## average across samples then take mean square across time steps
C = t(apply(C,1,prop))

## remove effects on bot
J[1,] = 0
C[1,] = 0

## thresholding
J = J*(C>0.1)
C = C*(C>0.1)

##
## DYNAMICAL INTERACTION PLOT (v1)
pdf(paste(pathToOut,"/fig_DIN_v1.pdf",sep=""),width=10,height=10)
.plot.DIN(J,C,colnames(TS)[-1])
dev.off()

##
## DYNAMICAL INTERACTION PLOT (v2)
pdf(paste(pathToOut,"/fig_DIN_v2.pdf",sep=""),width=12,height=12)
.plot.DIN2(J,C,colnames(TS)[-1])
dev.off()

## visualise interpolation
par(mfrow=c(N,N),mar=c(1,1,1,1)*1,oma=c(2,2,1,1))
color_vector = rainbow(N,alpha=0.75)
for(i in 1:N)
{
    for(j in 1:N)
    {
        if (i == j)
        {
            plot(-1:1, -1:1, col="white", cex = 0, xaxt="n", yaxt="n", type="l")
            text(0, 0, colnames(TS)[-1][i], cex = 2)
        } else
        {
            x = apply(Yhat_o[[i]],2,mean)
            y = apply(Yhat_o[[j]],2,mean)
            plot(x,y,col=color_vector[i],type="l")
            points(TS[,i+1],TS[,j+1],pch=16)
            lines(x,y,col=color_vector[i],type="l")    
        }
    }
}
par(mfrow=c(1,1))

#
###

####################################
## CROSS-VALIDATION PROCESS MODEL ##
####################################

## goal: perform cross-validation of process model by predicting the second half of the time series (of the interpolated p.c. growth rates)
##       the cross-validation is performed to identify the best value for regularisation on the nonlinear part of the model

## load model_o
load(file=paste(pathToOut,"/","Yhat_o.RData"    ,sep=""))
load(file=paste(pathToOut,"/","ddt.Yhat_o.RData",sep=""))
load(file=paste(pathToOut,"/","Omega_o.RData"   ,sep=""))

## parameters for cross-validation
K_p             = 3                                        # number of models to fit per folds and regularisation parameter
max_fold        = 2/3                                      # beginning of test set
folds           = list(c(0,1/3) * max_fold, 
                       c(1/3,2/3) * max_fold, 
                       c(2/3,1) * max_fold)  # proportion of the data that should be considered for training and validation
crossValParVect = seq(0.005,0.05,0.005)

## run cross-validation
resultsCrossVal_p = crossVal_p(TS,alpha_i,Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,folds,crossValParVect)

## visualise
pdf(paste(pathToOut,"/fig_crossVal_p.pdf",sep=""))
plotCrossVal_p(resultsCrossVal_p)
dev.off()

## store results
save(resultsCrossVal_p,file=paste(pathToOut,"/","crossVal_p.RData"   ,sep=""))

#
###
