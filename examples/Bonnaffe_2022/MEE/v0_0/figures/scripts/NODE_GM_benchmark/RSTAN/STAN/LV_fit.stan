
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

  
