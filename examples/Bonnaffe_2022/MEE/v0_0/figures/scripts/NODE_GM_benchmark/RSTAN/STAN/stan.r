############
## stan.r ##
############

## goal: bayesian inference of the ecological dynamics of the Lynx hare system

##
library(deSolve)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##############
## INITIATE ##
##############

##
Y <- read.table("TS_3DLV.csv",sep=";",header=T)
Y <- Y[20:50,]
Y$t <- Y$t - min(Y$t) + 1
# Y[,2:4] <- t(t(Y[,2:4])/apply(Y[,2:4],2,max)*3)


#
###

###############
## PLOT DATA ##
###############

##
par(mfrow=c(1,1))
#
plot(Y$t, Y$B, type="b", lty=2, xlab="Time", ylab="Count", bty="n", pch=16, col="red",ylim=c(0,max(Y[,-1])))
lines(Y$t, Y$G, type="b", lty=2, pch=16,col="green")
lines(Y$t, Y$B, type="b", lty=2, pch=16,col="blue")
lines(Y$t, Y$R, type="b", lty=2, pch=16, col="red")
#
par(mfrow=c(1,1))

#
###

####################
## MODEL AND DATA ##
####################

## initial conditions and parameters
params <- c(runif(9,-0.01,0.01))
inits <- c(G=Y[1,2], B=Y[1,3], R=Y[1,4])
t_min <- 0
t_max <- 30
times <- seq(t_min,t_max,1)

## lotka-volterra model
dydt <- function(t, y, params) {
  list(c(
    y[1] * (params[1] + params[2] * y[2] + params[3] * y[3]),
    y[2] * (params[4] + params[5] * y[1] + params[6] * y[3]),
    y[3] * (params[7] + params[8] * y[1] + params[9] * y[2])
  ))
}

# run
out <- ode(inits, times, dydt, params, method="ode45")

# plot the system
out <- data.frame(out)
plot(NA,NA, xlim = c(t_min, t_max), ylim=c(0, max(out[,-1])), xlab = "Time", ylab="Fraction of Host Population")
lines(out$G ~ out$time, col="green")
lines(out$B ~ out$time, col="blue")
lines(out$R ~ out$time, col="red")
legend("topright", legend = c("G", "B", "R"), col = c("green", "blue", "red"), lty = c(1, 1, 1), bty="n")

#
###


###############
## INFERENCE ##
###############

## stan model
cat(
  '
  // ODE SYSTEM
  functions {
  
  real[] LV(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  
  real dydt[3];
  dydt[1] = y[1]*(params[1] + params[2]* y[2] + params[3] * y[3]);
  dydt[2] = y[2]*(params[4] + params[5]* y[1] + params[6] * y[3]);
  dydt[3] = y[3]*(params[7] + params[8]* y[1] + params[9] * y[2]);
  return dydt;

  }
  
  }
  
  // DATA
  data {

  int<lower = 1> n_params;
  real t0; // initial time  
  int<lower = 1> T; // final time
  real ts[T]; // vector of time steps
  real y[T,3]; // data time series
  
  
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  }
  
  // PARAMETERS
  parameters {

  real<lower = 0> y0[3];  
  real params[n_params];

  }
  
  transformed parameters {
  
  // model predictions
  real y_hat[T, 3];
  y_hat = integrate_ode_rk45(LV, y0, t0, ts, params, x_r, x_i);

  }

  // MODEL
  model {

  // likelihood
  for(t in 1:T){
  y[t,1] ~ normal(log(y_hat[t,1]), 0.5);
  y[t,2] ~ normal(log(y_hat[t,2]), 0.5);
  y[t,3] ~ normal(log(y_hat[t,3]), 0.5);
  }

  // priors
  params ~ normal(0, 1000); //
  y0 ~ normal(0, 1000); //
  
  }

  ', 
  file = "LV_fit.stan", sep="", fill=T)

## data
stan_d = list(n_params = length(params),
              t0 = 0,
              T = max(Y$t),
              ts = Y$t,
              y = Y[,2:4])

## outputVect
outputVect = c("y_hat", "y0", "params")

##
initf1 <- function(n_chains) {
  inits <- list()
  for(i in 1:n_chains)
  {
    theta_0 = c(as.numeric(Y[1,-1]),runif(9,-0.01,0.01))
    names(theta_0) = c("G","B","R",
                       "beta_1","beta_2","beta_3",
                       "beta_4","beta_5","beta_6",
                       "beta_7","beta_8","beta_9")  
    inits[[i]] <- as.list(theta_0)
  }
  
  return(inits)
}

## Test / debug the model:
test = stan("LV_fit.stan",
            data = stan_d,
            pars = outputVect,
            init = initf1(1),
            chains = 1, 
            iter = 100)

## You should do some MCMC diagnostics, including:
traceplot(test, pars="lp__")
traceplot(test, pars=c("params", "y0"))
summary(test)$summary[,"Rhat"]

## Fit and sample from the posterior
mod = stan(fit = test,
           data = stan_d,
           pars = outputVect,
           init = initf1(1),
           chains = 1,
           warmup = 1000,
           iter = 10000)

# You should do some MCMC diagnostics, including:
traceplot(mod, pars="lp__")
traceplot(mod, pars=c("params", "y0"))
summary(mod)$summary[,"Rhat"]

# Extract the posterior samples to a structured list:
posts <- extract(mod)

#
###

############################
## post predictive checks ##
############################

##
names(posts)
hist(posts$params[,3])
plot(data.frame(cbind(posts$y0,posts$params))[seq(1,1000,5),])

##
MaP <- c(posts$y0[which.max(posts$lp__),], posts$params[which.max(posts$lp__),]); 
names(MaP) <- c("G","B","R",
                "beta_1","beta_2","beta_3",
                "beta_4","beta_5","beta_6",
                "beta_7","beta_8","beta_9")  

# run
out <- ode(MaP[1:3], times, dydt, MaP[-c(1:3)], method="ode45")

# plot the system
out <- data.frame(out)
plot(NA,NA, xlim = c(t_min, t_max), ylim=c(0, max(out[,-1])), xlab = "Time", ylab="Fraction of Host Population")
points(Y$t,Y$G,col="green")
points(Y$t,Y$B,col="blue")
points(Y$t,Y$R,col="red")
lines(out$G ~ out$time, col="green")
lines(out$B ~ out$time, col="blue")
lines(out$R ~ out$time, col="red")
legend("topright", legend = c("G", "B", "R"), col = c("green", "blue", "red"), lty = c(1, 1, 1), bty="n")

## write the chains results
chains <- cbind(posts$lp__,posts$y0,posts$params); colnames(chains) <- c("LogPost",names(MaP))
write.table(chains,"samples.csv", sep=";", row.names = FALSE)

#
###
