#################
## f_NODE_GM.r ##
#################

## goal: supporting functions to implement gradient matching (GM) fitting of neural ordinary differential equations (NODEs)

## author: Willem Bonnaffe (w.bonnaffe@gmail.com)

## method:
## 1. interpolate the time series with sin ANN functions
## 2. estimate linear and non-linear coupling between the time series

## update log:
## 09-06-2022 - created v0_0

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
