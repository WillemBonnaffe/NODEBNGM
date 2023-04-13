#########################
## f_NODE_GM_RStudio.r ##
#########################

## goal: supporting functions to implement gradient matching (GM) fitting of neural ordinary differential equations (NODEs)
##       the repository is designed with RStudio in mind 
##       see associated script for implementing the approach ("m0_main_RStudio.r")

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## method:
## 1. interpolate the time series with neural networks with sinusoid activation functions
## 2. estimate linear and non-linear coupling between the time series with neural networks

## update log:
## 09-06-2022 - created v0_0
## 31-10-2022 - created v0_1
##            - corrected indexing of figures
## 12-04-2023 - implemented more than two-fold validation

###################
## FUNCTIONS SLP ##
###################

## goal: functions to define a single layer perceptron

## activation functions ##
lin       = function(x) x
ddx.lin   = function(x) 1
pol2      = function(x) x^2
ddx.pol2  = function(x) 2*x
sigmo     = function(x) 1/(1+exp(-x)) 
ddx.sigmo = function(x) sigmo(x) * (1 - sigmo(x))
expo      = function(x) exp(x)
ddx.expo  = function(x) exp(x)
relu      = function(x) (x>0)*x
ddx.relu  = function(x) (x>0)
tanh      = function(x) (exp(2*x)-1)/(exp(2*x)+1) 
ddx.tanh  = function(x) (2*exp(2*x))/(exp(2*x)+1) + (exp(2*x)-1)*(-1)*(2*exp(2*x))/((exp(2*x)+1)^2)

## SLP ##
## goal: single layer perceptron returning single output with input x and parameter vector Omega 
# x       - vector - input variables
# Omega   - vector - parameters
# f_sigma - func - activation function
SLP = function(x,Omega,f_sigma) 
{    
    Omega = matrix(Omega,ncol=2 + length(x))
    return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
}

## ddOmega.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each parameter 
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddOmega.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
    Omega     = matrix(Omega,ncol=2 + length(x))
    x         = t(x)
    Omega_1   = Omega[,1]
    Omega_2   = Omega[,2]
    Omega_3   = Omega[,-c(1:2)]
    ddOmega_1 = f_sigma(Omega_2 + Omega_3%*%t(x))
    ddOmega_2 = Omega_1 * ddu.f_sigma(Omega_2 + Omega_3%*%t(x))
    ddOmega_3 = Omega_1%*%x * as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x)))
    return(c(ddOmega_1,ddOmega_2,ddOmega_3))
}

## ddx.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each input variable
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddx.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
    Omega = matrix(Omega,ncol=2 + length(x))
    x         = t(x)
    Omega_1   = Omega[,1]
    Omega_2   = Omega[,2]
    Omega_3   = Omega[,-c(1:2)]
    ddx = Omega_1%*%(Omega_3*as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x))))
    return(ddx)
}

#
###

##############################
## FUNCTIONS BAYESIAN MODEL ##
##############################

## goal: functions to define a simple Bayesian model with Gaussian error structure

## r2_p ##
## goal: compute r squarred of predictions given data
# X     - matrix - explanatory variables
# Y     - vector - response variable
# f     - func   - predictive function 
# Omega - vector - parameters of predictive function
r2 = function(X,Y,f,Omega)
{
  res      = Y - f(X,Omega) 
  r2       = 1 - sd(res^2)/sd(Y^2)
  return(r2)
}

## logLik ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
logLik = function(X,Y,f,Omega,sd_1)
{
  res      = Y - f(X,Omega)
  logLik   = - sum((res^2)/(sd_1^2))
  return(logLik)
}

## logPrior ##
## goal: compute log prior density of the process model 
# Omega - vector - parameters
# sd_2  - float  - standard deviation of prior
logPrior = function(Omega,sd_2)
{
  logPrior = - sum((Omega^2)/(sd_2^2))
  return(logPrior)
}

## logPost ##
## goal: log posterior distribution with normal error 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
logPost = function(X,Y,f,Omega,sd_1,sd_2)
{
    res      = Y - f(X,Omega)
    logLik   = - sum(  (res^2)/(sd_1^2))
    logPrior = - sum((Omega^2)/(sd_2^2))
    logPost  = logLik + logPrior
    return(logPost)
}

## ddOmega.logPost ##
## goal: compute the derivate of the log posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
ddOmega.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
    res              = Y - f(X,Omega)
    ddOmega.res      =   - df(X,Omega)
    ddOmega.logLik   = - 2 * ddOmega.res%*%res/(sd_1^2)
    ddOmega.logPrior = - 2 * Omega/(sd_2^2)
    ddOmega.logPost  = ddOmega.logLik + ddOmega.logPrior
    return(ddOmega.logPost)
}

## argmax.logPost ##
## goal: compute parameter vector that maximises log posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
argmax.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
    error_     = function(x) -logPost(X,Y,f,x,sd_1,sd_2)
    graderror_ = function(x) -ddOmega.logPost(X,Y,f,df,x,sd_1,sd_2)
    Omega      = optim(par    = Omega,
                       fn     = error_,
                       gr     = graderror_,
                       method = "BFGS"#,
                       # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
                       )$par
    return(Omega)
}

## logMarLik ##
## goal: compute the log marginal likelihood 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarLik = function(X,Y,f,Omega)
{
    res      = Y - f(X,Omega)
    logMarLik   = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
    return(logMarLik)
}

## logMarPri ##
## goal: compute the log marginal prior density
# Omega - vector - parameters
logMarPri = function(Omega)
{
    logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1) 
    return(logMarPri)
}

## logMarPost ##
## goal: compute the log marginal posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
logMarPost = function(X,Y,f,Omega,c=1)
{
    res       = Y - f(X,Omega)
    logMarLik = - 0.5 * length(Y)     * log(0.5 * sum(res^2)   + 1)
    logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
       logMarPos = logMarLik + c*logMarPri
    return(logMarPos)
}

## ddOmega.logMarPost ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
ddOmega.logMarPost = function(X,Y,f,df,Omega,c=1)
{
res                = Y - f(X,Omega)
    ddOmega.res        =   - df(X,Omega)
    ddOmega.logMarLik  = - 0.5 * length(Y)     * 1/(0.5 * sum(res^2)   + 1) * 0.5 * ddOmega.res%*%res
    ddOmega.logMarPri  = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
       ddOmega.logMarPos  = ddOmega.logMarLik + c*ddOmega.logMarPri ## divide by number of neurones in the network
    return(ddOmega.logMarPos)
}

## argmax.logMarPost ##
## goal: compute parameter vector that maximises the log marginal density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# c     - scalar - parameter to control strength of regularisation
argmax.logMarPost = function(X,Y,f,df,Omega,c=1)
{
    error_     = function(x) -logMarPost(X,Y,f,x,c)
    graderror_ = function(x) -ddOmega.logMarPost(X,Y,f,df,x,c)
    Omega      = optim(par    = Omega,
                       fn     = error_,
                       gr     = graderror_,
                       method = "BFGS"# ,
                       # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
                       )$par
    return(Omega)
}

#
###

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## loadData_o ##
## goal: prepare data for training of observation model
# TS      - matrix - time series (time in first column, state variables in remaining columns)
# alpha_i - scalar - interpolation factor (e.g. 2 means that each time step is divided by 2)
loadData_o = function(TS,alpha_i)
{
    ## data specs
    N       = ncol(TS) - 1
    n       = nrow(TS)
        
    ## predictive and response variable
    t   = TS[,1] 
    Y   = TS[,-1]
    nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)
    
    ## standardise time steps
    t_0 = min(t)
    t_f = max(t)
    dt  = diff(t[1:2])
    t_  = 1:length(t)
    nt_ = seq(min(t_),max(t_),(t_[2]-t_[1])/alpha_i)
    
    ## standardise data
    Y_     = log(Y)
    mean_y = apply(Y_,2,mean)
    sd_y   = apply(Y_,2,sd)
    Y_     = t((t(Y_)-mean_y)/sd_y)

    ## terminate
    return(list("t"      = t,
                "t_"     = t_,
                "nt"     = nt,
                "nt_"    = nt_,
                "Y"      = Y,
                "Y_"     = Y_,
                "mean_y" = mean_y,
                "sd_y"   = sd_y,
                "dt"     = dt,
                "N"      = N
                ))
}

## f_o ##
## goal: compute predicted values of response variable at time step t
# t     - float  - time step 
# Omega - vector - parameters 
f_o = function(t,Omega)
{    
    Omega = matrix(Omega,ncol=3)
    return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o ##
## goal: compute time derivative of the predicted response t time step t
# t     - float  - time step 
# Omega - vector - parameters 
ddt.f_o = function(t,Omega)
{    
    Omega = matrix(Omega,ncol=3)
    return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega.f_o ##
## goal: compute derivative of the predicted response wtr to each network parameter
# t     - float  - time step
# Omega - vector - parameters 
ddOmega.f_o = function(t,Omega)
{    
    Omega      = matrix(Omega,ncol=3)
    dfdOmega_1 =                sin(pi*(t*Omega[,2] + Omega[,3]))
    dfdOmega_2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
    dfdOmega_3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
    return(c(dfdOmega_1,dfdOmega_2,dfdOmega_3))
}

## *.eval ##
## goal: compute functions across multiple time steps
# t     - vector - time steps in arbitrary units
# Omega - vector - parameters 
f_o.eval         = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval     = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega.f_o(x,Omega))

## trainModel_o ##
## goal: train observation model
# TS      - matrix - time series (time in first column, state variables in remaining columns)
# alpha_i - int    - interpolation factor (e.g. 2 means that each time step is divided by 2)
# N_o     - int    - number of parameters in the observation model
# K_o     - int    - number of samples to take
# rho     - float  - [0,1] proportion of best samples to keep or reject
trainModel_o = function(TS,alpha_i,N_o,K_o,rho=1)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)

    ## interpolate each column in the time series
    Omega_o    = list()
    Yhat_o     = list()
    ddt.Yhat_o = list()
    logPost_o  = list()
    for (i in 1:N)
    {
         ## iterator
         message(paste("fitting: ",i,"/",N,sep=""))
    
        ## get ensemble
        Omega_o_i = NULL
        logPost_o_i = NULL
        for(k in 1:K_o)
        {   
            check = F
            while(check == F)
            {
                ## fit
                Omega_0   = rnorm(N_o[i],0,0.001)
                Omega_f   = argmax.logMarPost(t_,Y_[,i],f_o.eval,ddOmega.f_o.eval,Omega_0,3/N_o[i])
    
                ## update
                logPost_0 = logMarPost(t_,Y_[,i],f_o.eval,Omega_0,3/N_o[i])
                logPost_f = logMarPost(t_,Y_[,i],f_o.eval,Omega_f,3/N_o[i])
                message(paste(k,"/",K_o,"\t",
                        format(round(logPost_0,2),nsmall=2),"\t","-->","\t",
                        format(round(logPost_f,2),nsmall=2),sep=""))
                
                ## reject or keep sample
                check = (logPost_f >= logPost_0 + 10)
                if(check == T)
                {
                       Omega_o_i = rbind(Omega_o_i,Omega_f)
                    logPost_o_i = c(logPost_o_i,logMarLik(t_,Y_[,i],f_o.eval,Omega_f))
                }
            }
        }
    
        ## sort by best models
        s = rev(order(logPost_o_i))[1:round(rho*K_o)]
        Omega_o_i = Omega_o_i[s,]
    
        ## store ensemble
        Omega_o[[i]] = Omega_o_i
    
        ## compute predictions
        Yhat_o[[i]]     = t(apply(Omega_o[[i]],1,function(x)     f_o.eval(nt_,x))) 
        ddt.Yhat_o[[i]] = t(apply(Omega_o[[i]],1,function(x) ddt.f_o.eval(nt_,x)))
    
        ## de-standardise
        Yhat_o[[i]]     = exp(mean_y[i] + sd_y[i] * Yhat_o[[i]])
        ddt.Yhat_o[[i]] = 1/dt * sd_y[i] * Yhat_o[[i]] * ddt.Yhat_o[[i]]
     }
    
    ## store results
    names(Yhat_o)     = colnames(TS[,-1])
    names(ddt.Yhat_o) = colnames(TS[,-1])
    return(list("Yhat_o"     = Yhat_o,
                "ddt.Yhat_o" = ddt.Yhat_o,
                "Omega_o"    = Omega_o))
}

## plot_o ## 
## goal: visualise observation model fit
# TS          - matrix - time series (time in first column, state variables in remaining columns)
# alpha_i     - int    - interpolation factor (e.g. 2 means that each time step is divided by 2)
# Yhat_o      - list   - list containing interpolation ensembles of each state variables
# ddt.Yhat_o  - list   - list containing temporal derivative of interpolation ensembles of each state variables
plotModel_o = function(TS,alpha_i,Yhat_o,ddt.Yhat_o)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)

    col   = (rainbow(N,start=0.1,end=0.9))
    xlab  = c("","Time")
    ylab  = c("Y(t)","dY/dt(t)")
    index = NULL; for (i in 1:N) index = c(index,paste(i,c("a.","b."),sep="")) 
    par(mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
    layout(cbind(c(1,2),c(3,4),c(5,6))[,1:min(3,N)])
    for (i in 1:N)
    {
        ## predictions 
        E.Yhat_o       = apply(Yhat_o[[i]],2,mean)
        q05.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.05)
        q95.Yhat_o     = apply(Yhat_o[[i]],2,quantile,p=0.95)
        E.ddt.Yhat_o   = apply(ddt.Yhat_o[[i]],2,mean)
        q05.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.05)
        q95.ddt.Yhat_o = apply(ddt.Yhat_o[[i]],2,quantile,p=0.95)
        #
        ## training data
        t = TS[,1] 
        Y = TS[,-1][,i] 
        nt = seq(min(t),max(t),(t[2]-t[1])/alpha_i)
        #
        ## visualise interpolated response
        plot(t,rep(0,length(t)),ylim=c(min(Y),max(Y))+c(-1,1)*0.2*(max(Y)-min(Y)),type="l",lty=3,xlab=xlab[1],ylab=if(i%%3 == 1)ylab[1]else"")
        polygon(c(nt,rev(nt)),c(q05.Yhat_o,rev(q95.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
        points(t,Y,pch=16,col=adjustcolor("black",0.75)) 
        lines(nt,E.Yhat_o,col=adjustcolor(col[i],0.75),lwd=2)
        if(!is.null(index)) legend("topright",legend=index[1+(i-1)*2],bty="n",cex=1.5)
        legend("bottom",legend=colnames(TS)[-1][i],lty=1,col=col[i],lwd=2,bty="n",horiz=T)
        #
        ## visualise temporal derivative
        plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddt.Yhat_o))*2,type="l",lty=3,xlab=xlab[2],ylab=if(i%%3 == 1)ylab[2]else"")
        polygon(c(nt,rev(nt)),c(q05.ddt.Yhat_o,rev(q95.ddt.Yhat_o)),col=adjustcolor(col[i],alpha=0.25),border=NA)
        lines(nt,E.ddt.Yhat_o,col=adjustcolor(col[i],0.75),lwd=2)
        if(!is.null(index)) legend("topright",legend=index[2+(i-1)*2],bty="n",cex=1.5)
        legend("bottom",legend=colnames(TS)[-1][i],lty=1,col=col[i],lwd=2,bty="n",horiz=T)
    }
    par(mfrow=c(1,1))
}

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## loadData_p ##
## goal: prepare data for training process model
# Yhat_o      - list   - list containing interpolation ensembles of each state variables
# ddt.Yhat_o  - list   - list containing temporal derivative of interpolation ensembles of each state variables
loadData_p = function(Yhat_o,ddt.Yhat_o)
{
    ## data specs
    N   = length(Yhat_o)
    n   = ncol(Yhat_o[[1]])
    K_o = nrow(Yhat_o[[1]])

    ## predictive variable
    E.X_p     = NULL
    E.ddt.X_p = NULL
    E.Y_p     = NULL
    X_p       = NULL
    ddt.X_p   = NULL
    Y_p       = NULL
    for(j in 1:N)
    {
        # ## select random ensemble element (for stochastic training)
        # v       = sample(1:K_o,1)
        # X_p     = cbind(X_p,        Yhat_o[[j]][v,])
        # ddt.X_p = cbind(ddt.X_p,ddt.Yhat_o[[j]][v,])
        # Y_p     = cbind(Y_p,1/Yhat_o[[j]][v,]*ddt.Yhat_o[[j]][v,])
    
        ## mean of ensemble
        E.X_p     = cbind(E.X_p,     apply(Yhat_o[[j]],2,mean))
        E.ddt.X_p = cbind(E.ddt.X_p, apply(ddt.Yhat_o[[j]],2,mean))
        E.Y_p     = cbind(E.Y_p,     apply(1/Yhat_o[[j]]*ddt.Yhat_o[[j]],2,mean))
    }
    
    ## variables
    X     = E.X_p
    ddt.X = E.ddt.X_p
    Y     = E.Y_p
    
    ## standardise predictive variables
    X_     = X
    mean_x = apply(E.X_p,2,mean)
    sd_x   = apply(E.X_p,2,sd)
    X_     = t((t(X_)-mean_x)/sd_x)
    
    ## standardise response variable
    Y_     = Y
    mean_y = apply(E.Y_p,2,mean)
    sd_y   = apply(E.Y_p,2,sd)
    Y_     = t((t(Y_)-mean_y)/sd_y)

    ## terminate
    return(list("X"      = X,
                "ddt.X"  = ddt.X,
                "X_"     = X_,
                "Y"      = Y,
                "Y_"     = Y_,
                "mean_x" = mean_x,
                "sd_x"   = sd_x,
                "mean_y" = mean_y,
                "sd_y"   = sd_y
                ))
}

## activation function of process model ## 
# f_sigma_p     = lin
# ddx.f_sigma_p = ddx.lin
# f_sigma_p     = pol2   
# ddx.f_sigma_p = ddx.pol2
f_sigma_p     = expo
ddx.f_sigma_p = ddx.expo
# f_sigma_p     = sigmo     
# ddx.f_sigma_p = ddx.sigmo 

## f_p ##
## goal: compute predicted response variable of process model defined given vector of explanatory variable at a given time step 
# x       - vector - input variables
# Omega   - vector - parameters
f_p = function(x,Omega)
{
    Omega = matrix(Omega,ncol=2)
    return(SLP(x,Omega[,1],lin) + SLP(x,Omega[,2],f_sigma_p))
}

## ddOmega.f_p ##
## goal: compute derivative vector of the process model wtr to each parameter at a given time step
# x      - vector - input variables
# Omega  - vector - parameters 
ddOmega.f_p = function(x,Omega)
{
    Omega = matrix(Omega,ncol=2)
    return(c(ddOmega.SLP(x,Omega[,1],lin,ddx.lin),ddOmega.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p)))
}

## ddx.f_p ##
## goal: compute derivative vector of the process model wtr to each state variable at a given time step
# x      - vector - input variables
# Omega  - vector - parameters
ddx.f_p = function(x,Omega)
{
    Omega = matrix(Omega,ncol=2)
    return(ddx.SLP(x,Omega[,1],lin,ddx.lin) + ddx.SLP(x,Omega[,2],f_sigma_p,ddx.f_sigma_p))
}

## *.eval ##
## goal: compute process model functions across several time steps
# X      - matrix - matrix of variables across all time steps
# ddt.X  - matrix - matrix of temporal derivatives of variables
# Omega  - vector - vector of parameters of the network 
f_p.eval           = function(X,Omega) apply(t(X),2,function(x) f_p(x,Omega))
ddx.f_p.eval       = function(X,Omega) apply(t(X),2,function(x) ddx.f_p(x,Omega))
ddOmega.f_p.eval   = function(X,Omega) apply(t(X),2,function(x) ddOmega.f_p(x,Omega))

## trainModel_p ##
## goal: train process model
# Yhat_o      - list   - list containing interpolation ensembles of each state variables
# ddt.Yhat_o  - list   - list containing temporal derivative of interpolation ensembles of each state variables
# N_p         - int    - number of parameters in process model
# sd1_p       - float  - standard deviation of the likelihood
# sd2_p       - vector - standard deviation of the prior distributions
# K_p         - int    - number of elements to sample
# trainSplit  - float  - [0,1] proportion of the data to use for training vs testing
trainModel_p = function(Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,trainSplit=0.75)
{
    ## for each variable
    Omega_p    = list()
    Yhat_p     = list()
    ddx.Yhat_p = list()
    Geber_p    = list()
    for(i in 1:N)
    {
        ## iterator
        message(paste("fitting: ",i,"/",N,sep=""))
    
        ## fit model
        Omega_p_i = NULL
        for(k in 1:K_p)
        {
    
            ## dataloader
            attach(loadData_p(Yhat_o,ddt.Yhat_o),warn.conflicts=F)
 
            ## split train and test
            s_l = 1:round(nrow(X_)*trainSplit)
            s_t = -s_l
            
            ## fit
            Omega_0      = rnorm(N_p[i],0,0.001)
            Yhat         = function(X,Omega) f_p.eval(X,Omega)
            ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega)
            Omega_f      = argmax.logPost(X_[s_l,],Y_[s_l,i],Yhat,ddOmega.Yhat,Omega_0,sd1_p,sd2_p[[i]])
            Omega_p_i    = rbind(Omega_p_i,Omega_f)
    
            ## update
            logPost_0_l    = logLik(X_[s_l,],Y_[s_l,i],Yhat,Omega_0,sd1_p)
            logPost_f_l    = logLik(X_[s_l,],Y_[s_l,i],Yhat,Omega_f,sd1_p)
            logPost_0_t    = logLik(X_[s_t,],Y_[s_t,i],Yhat,Omega_0,sd1_p)
            logPost_f_t    = logLik(X_[s_t,],Y_[s_t,i],Yhat,Omega_f,sd1_p)
            message(paste(k,"/",K_p,"\t",
                    "initial_train: ",format(round(logPost_0_l,2),nsmall=2),"\t",
                    "final_train: ",format(round(logPost_f_l,2),nsmall=2),"\t",
                    "initial_test: ",format(round(logPost_0_t,2),nsmall=2),"\t",
                    "final_test: ",format(round(logPost_f_t,2),nsmall=2),sep=""))
        }
    
        ## store
        Omega_p[[i]] = Omega_p_i
    
        ## predictions
        Yhat_p[[i]]     = t(apply(Omega_p[[i]],1,function(x) mean_y[i] + sd_y[i] * f_p.eval(X_,x)))
        ddx.Yhat_p[[i]] = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * ddx.f_p.eval(X_,x)))
        Geber_p[[i]]    = t(apply(Omega_p[[i]],1,function(x) sd_y[i] * 1/sd_x * t(ddt.X) * ddx.f_p.eval(X_,x)))
    }

    ## terminate
    return(list("Yhat_p"     = Yhat_p,
                "ddx.Yhat_p" = ddx.Yhat_p,
                "Geber_p"    = Geber_p,
                "Omega_p"    = Omega_p
                ))
}

## plotModel_p ##
# TS          - matrix - time series (time in first column, state variables in remaining columns)
# alpha_i     - int    - interpolation factor (e.g. 2 means that each time step is divided by 2)
# Yhat_p      - list   - list (for each state variable) containing ensembles of predicted response by the process model 
# ddx.Yhat_p  - list   - list (for each state variable) containing ensembles of derivative of predicted response wtr to inputs
# Geber_p     - list   - list (for each state variable) containing ensembles of contributions (see Geber method)
plotModel_p = function(TS,alpha_i,Yhat_p,ddx.Yhat_p,Geber_p)
{
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)
    attach(loadData_p(Yhat_o,ddt.Yhat_o),warn.conflicts=F)
    
    col   = (rainbow(N,start=0.1,end=0.9))
    xlab  = c("","","Time")
    ylab  = c("P.c. growth rate","Effects","Contributions")
    index = NULL; for (i in 1:N) index = c(index,paste(i,c("a.","b.","c."),sep="")) 
    legend = paste(colnames(TS)[-1])
    par(mar=c(4,4.5,0,0),oma=c(1,1,1,3),cex.lab=1.5)
    layout(cbind(c(1,2,3),c(4,5,6),c(7,8,9))[,1:min(3,N)])
    for(i in 1:N)
    {
        ## predictions
        E.Yhat_p       = apply(Yhat_p[[i]],2,mean)
        q05.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.05)
        q95.Yhat_p     = apply(Yhat_p[[i]],2,quantile,p=0.95)
        E.ddx.Yhat_p   = t(matrix(apply(ddx.Yhat_p[[i]],2,mean),ncol=nrow(X)))
        q05.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
        q95.ddx.Yhat_p = t(matrix(apply(ddx.Yhat_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
        E.Geber_p      = t(matrix(apply(Geber_p[[i]],2,mean),ncol=nrow(X)))
        q05.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.05),ncol=nrow(X)))
        q95.Geber_p    = t(matrix(apply(Geber_p[[i]],2,quantile,p=0.95),ncol=nrow(X)))
        #
        ## dynamics
        x = nt
        y = Y[,i]
        plot(x,rep(0,length(x)),ylim=c(-1,1)*max(abs(y))*1.5,type="l",lty=3,xlab=xlab[1],ylab=if(i%%3==1)ylab[1]else"")
        points(x,y,pch=16,col=adjustcolor("black",0.75))
        polygon(c(nt,rev(nt)),c(q05.Yhat_p,rev(q95.Yhat_p)),col=adjustcolor(col[i],alpha=0.2),border=NA)
        lines(nt,E.Yhat_p,col=adjustcolor(col[i],alpha=0.75),lwd=2)
        if(!is.null(index)) legend("topright",legend=index[1+(i-1)*(3)],bty="n",cex=1.5)
        #
        ##
        if(!is.null(legend)) legend("bottom" ,legend=legend[i],pch=15,col=adjustcolor(col,alpha=0.75)[i],bty="n",horiz=F,cex=0.75)
        #
        ## effects
        plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.ddx.Yhat_p))*1.5,type="l",lty=3,xlab=xlab[2],ylab=if(i%%3==1)ylab[2]else"")
        for(j in 1:N) lines(nt,E.ddx.Yhat_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
        for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.ddx.Yhat_p[,j],rev(q95.ddx.Yhat_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
        if(!is.null(index))  legend("topright",legend=index[2+(i-1)*(3)],bty="n",cex=1.5)
        # if(!is.null(legend)) legend("bottom" ,legend=legend,lty=1,col=adjustcolor(col,alpha=0.75),bty="n",horiz=T,lwd=2)
        #
        ## legend
        if(!is.null(legend)) legend("bottomleft" ,legend=legend[1:(N/3)],pch=15,col=adjustcolor(col,alpha=0.75)[1:(N/3)],bty="n",horiz=F,cex=0.75)
        if(!is.null(legend)) legend("bottom"  ,legend=legend[((N/3)+1):(2*N/3)],pch=15,col=adjustcolor(col,alpha=0.75)[((N/3)+1):(2*N/3)],bty="n",horiz=F,cex=0.75)
        if(!is.null(legend)) legend("bottomright"  ,legend=legend[((2*N/3)+1):N],pch=15,col=adjustcolor(col,alpha=0.75)[((2*N/3)+1):N],bty="n",horiz=F,cex=0.75)
        #
        ## Geber
        plot(nt,rep(0,length(nt)),ylim=c(-1,1)*max(abs(E.Geber_p))*1.5,type="l",lty=3,xlab=xlab[3],ylab=if(i%%3==1)ylab[3]else"")
        for(j in 1:N) lines(nt, E.Geber_p[,j],col=adjustcolor(col[j],alpha=0.75),lwd=2)
        for(j in 1:N) polygon(c(nt,rev(nt)),c(q05.Geber_p[,j],rev(q95.Geber_p[,j])),col=adjustcolor(col[j],alpha=0.2),border=NA)
        if(!is.null(index))  legend("topright",legend=index[3+(i-1)*(3)],bty="n",cex=1.5)
        #
        ## legend
        if(!is.null(legend)) legend("bottomleft" ,legend=legend[1:(N/3)],pch=15,col=adjustcolor(col,alpha=0.75)[1:(N/3)],bty="n",horiz=F,cex=0.75)
        if(!is.null(legend)) legend("bottom"  ,legend=legend[((N/3)+1):(2*N/3)],pch=15,col=adjustcolor(col,alpha=0.75)[((N/3)+1):(2*N/3)],bty="n",horiz=F,cex=0.75)
        if(!is.null(legend)) legend("bottomright"  ,legend=legend[((2*N/3)+1):N],pch=15,col=adjustcolor(col,alpha=0.75)[((2*N/3)+1):N],bty="n",horiz=F,cex=0.75)
     }
    par(mfrow=c(1,1))
}

## crossVal_p ##
## goal: preform cross validation of process model by progressively decreasing the constraint on nonlinear component of model
# TS              - matrix - time series (time in first column, state variables in remaining columns)
# alpha_i         - int    - interpolation factor (e.g. 2 means that each time step is divided by 2)
# Yhat_o          - list   - list containing interpolation ensembles of each state variables
# ddt.Yhat_o      - list   - list containing temporal derivative of interpolation ensembles of each state variables
# N_p             - int    - number of parameters in process model
# sd1_p           - float  - standard deviation of the likelihood
# sd2_p           - vector - standard deviation of the prior distributions
# K_p             - int    - number of elements to sample
# folds           - list   - list containing vectors indicating different ways of splitting the data for training and testing
#                            for instance c(1/2,1) means testing on the second half of the data
# crossValParVect - vector - vector of values for the regularisation 
crossVal_p = function(TS,alpha_i,Yhat_o,ddt.Yhat_o,N_p,sd1_p,sd2_p,K_p,folds,crossValParVect)
{
  
    ## load data
    attach(loadData_o(TS,alpha_i),warn.conflicts=F)
    attach(loadData_p(Yhat_o,ddt.Yhat_o),warn.conflicts=F)
  
    ## data specs
    N   = length(Yhat_o)
    n   = ncol(Yhat_o[[1]])

    ## for each variable
    Omega_p    = list()
    crossVal_p = list()
    for(i in 1:N)
    {
        ## iterator
        message(paste("fitting: ",i,"/",N,sep=""))
    
        ## cross validation
        crossVal_i   = NULL
        Omega_p[[i]] = list()
        for(k in 1:length(crossValParVect)) # for each regularisation param
        {
            ## iterator
            message(paste("crossval: ",k,"/",length(crossValParVect),sep=""))
    
            ## multiple folds
            crossVal_ik          = NULL
            Omega_p[[i]][[k]]    = list()
            for(u in 1:length(folds)) # for each fold in the data
            {
                ## iterator
                message(paste("fold: ",u,"/",length(folds),sep=""))
    
                ## fit model
                Omega_p_iku = NULL
                for(m in 1:K_p)
                {
                    ## dataloader
                    attach(loadData_o(TS,alpha_i),warn.conflicts=F)
    
                    ## split training/test (v0_1)
                    s_val  = round(folds[[u]][1]*n+1):round(folds[[u]][2]*n)
                    X_l    = X_[-s_val,]
                    Y_l    = Y_[-s_val,]
                    X_t    = X_[s_val,]
                    Y_t    = Y_[s_val,]
    
                    ## TO MODULARISE ##
                    sd2_p[[i]] = crossValParVect[k] # for linear and nonlinear part of network
                    # sd2_p[[i]][(N_p[i]/2):N_p[i]] = crossValParVect[k] # for nonlinear part of network only
    
                    ## fit
                    Omega_0      = rnorm(N_p[i],0,0.001)
                    Yhat         = function(X,Omega) f_p.eval(X,Omega)
                    ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X,Omega)
                    # Omega_f      = argmax.logMarPost(X_l,Y_l[,i],Yhat,ddOmega.Yhat,Omega_0,1/W_p[i])
                    Omega_f      = argmax.logPost(X_l,Y_l[,i],Yhat,ddOmega.Yhat,Omega_0,sd1_p,sd2_p[[i]])
                    Omega_p_iku  = rbind(Omega_p_iku,Omega_f)
    
                    ## update
                    logPost_0 = logPost(X_l,Y_l[,i],Yhat,Omega_0,sd1_p,sd2_p[[i]])
                    logPost_f = logPost(X_l,Y_l[,i],Yhat,Omega_f,sd1_p,sd2_p[[i]])
    
                    ## model performance
                    logLik_l = logLik(X_l,Y_l[,i],Yhat,Omega_f,sd1_p)
                    logLik_t = logLik(X_t,Y_t[,i],Yhat,Omega_f,sd1_p)

                    ## cross validation matrix
                    crossVal_ik  = rbind(crossVal_ik,cbind(logLik_l,logLik_t))
                }
                Omega_p[[i]][[k]][[u]] = Omega_p_iku
            }
    
            ## store
            E.crossVal_ik  = apply(crossVal_ik,2,mean)
            sd.crossVal_ik = apply(crossVal_ik,2,sd)
            crossVal_i     = rbind(crossVal_i,c(crossValParVect[k],E.crossVal_ik,sd.crossVal_ik))
            message(paste("logLik l vs t: ",
                    format(round(E.crossVal_ik[1],2),nsmall=2),"\t",
                    format(round(E.crossVal_ik[2],2),nsmall=2),"\n",sep=""))
        }
        message("\n")
    
        ## store
        crossVal_p[[i]] = crossVal_i
        colnames(crossVal_p[[i]]) = c("sd","logLik_l","logLik_t","sd.logLik_l","sd.logLik_t")
        # colnames(crossVal_p[[i]]) = c("w","logMarLik_l","logMarLik_t","sd.logMarLik_l","sd.logMarLik_t")
    }
    message("\n")
    return(crossVal_p)
}

## figure cross validation ##
## goal: visualise the results of the cross validation
# resultsCrossVal - list - list containing results of the cross validation
plotCrossVal_p = function(resultsCrossVal)
{
    index = c("a.","b.","c.","d.","e.","f.","g.","h.","i.","j.","k.","l.","m.","n.","o.","p.","q.","r.","s.","t.","u.","v.","w.","x.","y.","z.")
    colVect = rainbow(2,start=0.1,end=0.9)
    par(mfrow=c(min(3,N),1),mar=c(4,4.5,0,0),oma=c(1,1,1,1),cex.lab=1.5)
    for(i in 1:N)
    {
        ## unpack
        crossVal = resultsCrossVal_p[[i]]
        #
        ## plot
        x    = crossVal[,"sd"]
        y    = crossVal[,c("logLik_l","logLik_t")]
        sd.y = crossVal[,c("sd.logLik_l","sd.logLik_t")]
        xlab = if(i == N) "Std. prior" else ""
        ylab = "Log likelihood"
        plot(x,rep(0,length(x)),ylim=c(min(y-sd.y)-0.5*abs(min(y-sd.y)),max(y+sd.y)+0.5*abs(max(y+sd.y))),type="l",lty=3,xlab=xlab,ylab=ylab)
        #
        ## training line
        col  = colVect[1]
        y    = crossVal[,"logLik_l"]
        sd.y = crossVal[,"sd.logLik_l"]
        lines(x,y,col=col)
        polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
        #
        ## testing line
        col  = colVect[2]
        y    = crossVal[,"logLik_t"]
        sd.y = crossVal[,"sd.logLik_t"]
        lines(x,y,col=col)
        polygon(x=c(x,rev(x)),y=c(y-sd.y,rev(y+sd.y)),col=adjustcolor(col,0.25),border=NA)
        #
        ## legend
        legend("bottom",legend=c("Training","Validation"),col=colVect,lty=1,bty="n",cex=1.5,horiz=T)
        legend("topright",legend=index[i],bty="n",cex=1.5)
    }
    par(mfrow=c(1,1))
}

#
###

####################
## VISUALISATIONS ##
####################

## .plot.DIN
## goal: plot the dynamical interaction network of the system
# effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
# weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
# labels     - vector - vector of the names of the variables in the matrix
plotDIN = function(effectsMat,weightsMat,labels)
{
  ## dimensions
  N = dim(effectsMat)[1]
  
  ## scale effects and contributions
  effectsMat = (effectsMat>0)*1
  # weightsMat = weightsMat/sum(weightsMat) # proportion of total change
  
  ## angles
  theta = seq(0,2*pi,(2*pi)/N)
  x = cos(theta)
  y = sin(theta)
  x_ = cos(theta+(2*pi/N*0.25))
  y_ = sin(theta+(2*pi/N*0.25))
  x__ = 1.25*cos(theta+(2*pi/N*0.125))
  y__ = 1.25*sin(theta+(2*pi/N*0.125))
  
  ## plot interactions
  plot(x=c(-1:1)*1.5,y=c(-1:1)*1.5,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
  color_vector = adjustcolor(c("green","red"),alpha=0.5)
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      color_ = if(effectsMat[i,j]>0){color_vector[1]}else{color_vector[2]}
      # points(x__[i],y__[i],cex=30/N)
      text(x__[i],y__[i],labels=labels[i])
      if(weightsMat[i,j]*10>0)
      {
        arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
      }
    }
  }
  points(x,y,pch=16)
  points(x_,y_,pch=4)
  # legend("bottomleft",legend=c("Ingoing effects","Outgoing effects"),pch=c(4,16),lty=c(NA,NA),col=c("black","black"),cex=1.25,bty="n",horiz=F,lwd=2)
  # legend("bottomright",legend=c("Positive effects","Negative effects"),pch=c(NA,NA),lty=c(1,1),col=c(color_vector[1],color_vector[2]),cex=1.25,bty="n",horiz=F,lwd=2)
  legend("bottom",legend=c("Ingoing effects","Outgoing effects","Positive effects","Negative effects"),pch=c(4,16,15,15),col=c("black","black",color_vector[1],color_vector[2]),cex=1.25,bty="n",horiz=T)
}

## .plot.DIN2
## goal: plot the dynamical interaction network of the system
# effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
# weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
# labels     - vector - vector of the names of the variables in the matrix
plotDIN2 = function(J,C,labels)
{
    ## graphical parameters
    xlab = "Mean total effects"
    ylab =  "Log mean total contributions"
    x = apply(J,2,mean)
    y = log(apply(C,2,mean))
    delta_x = abs((max(x)-min(x)))*0.1
    delta_y = abs((max(y)-min(y)))*0.1
    color_vector = adjustcolor(c("green","red"),alpha=0.5)
    #
    ## main plotting area
    plot(x,y,xlim=c(min(x)-delta_x,max(x)+delta_x*1.5), ylim=c(min(y)-delta_y,max(y)+delta_y*1.5),cex=0,bty="n",xlab=xlab,ylab=ylab,cex.lab=1.5)
    lines(c(0,0),c(-10,10),lty=2)
    #
    ## interactions
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            colour = if (J[i,j] >= 0) color_vector[1] else color_vector[2]
            line_width = C[i,j]*10
            x_i = x[i]
            x_j = x[j]
            y_i = y[i]
            y_j = y[j]
            # if(C[i,j]*10>0) lines(x=c(x_i-delta_x/2,x_j),y=c(y_i,y_j),col=colour,lwd=line_width)
            if(C[i,j]*10>0) arrows(x0=x_j, x1=x_i-delta_x/2, y0=y_j, y1=y_i, col=colour, lwd=line_width)
        }
    }
    #
    ## legend
    # points(x-delta_x/2,y,col=rainbow(N,start=0.1,end=0.9)[-1],pch=17,cex=2)
    # points(x,        y,col=rainbow(N,start=0.1,end=0.9)[-1],pch=16,cex=2)
    points(x-delta_x/2,y,pch=4,cex=2)
    points(x,        y,pch=16,cex=2)
    text(x+delta_x, y,labels = labels,cex=1.25)
    legend("bottom",horiz=T,legend=c("Ingoing effects","Outgoing effects","Positive effects","Negative effects"),pch=c(4,16,NA,NA),lty=c(NA,NA,1,1),col=c("black","black",color_vector[1],color_vector[2]),cex=1.5,bg="white",box.col="white")
}

## plot_time_series
## goal: plot time series data
## TS - matrix - time series data
plotTimeSeries = function(TS)
{
    ## initiate
    N = ncol(TS) - 1
    colvect = rainbow(N, start=0.1, end=0.9)
    par(mfrow=c(N + 1, N), mar=c(4,4,2,2))
    #
    ## plots
    for (i in 1:N) plot(TS[,1], TS[,i+1], xlab="t", ylab="Y(t)", col=colvect[i], pch=16)
    for (i in 1:N)
    {
        for (j in 1:N)
        {
            if (i == j) 
            {
                plot(-1:1, -1:1, cex=0, xlab="", ylab="", xaxt="n", yaxt="n") 
                text(0, 0, colnames(TS)[i+1], cex = 2)
            } else 
            {
                plot(TS[,i+1], TS[,j+1], type="l", xlab="", ylab="", bty="n")    
            }
        }
    }
    #
    ## terminate
    par(mfrow=c(1,1), mar=c(5,5,3,3))
}

## plot_interpolations
## goal: alternative visualisation of interpolations 
## TS     - matrix - time series dataset
## Yhat_o - matrix - interpolated time series
plotInterpolations = function(TS, Yhat_o)
{
    ## initiate
    N = ncol(TS) - 1
    par(mfrow=c(N,N), mar=c(1,1,1,1)*1, oma=c(2,2,1,1))
    color_vector = rainbow(N, alpha=0.75)
    #
    ## plots
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
                lines(TS[,i+1],TS[,j+1],pch=16,type="b")
                lines(x,y,col=color_vector[i],type="l")    
            }
        }
    }    
    #
    ## terminate
    par(mfrow=c(1,1))
}

#
###