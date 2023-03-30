data {
  // This block declares all data which will be passed to the Stan model.
  // Pre period
  int<lower=0> N_pre;        // number of data items
  vector[N_pre] x_pre;      // ambient temperature
  vector[N_pre] y_pre;      // outcome energy vector
  // Post period
  int<lower=0> N_post;        // number of data items
  vector[N_post] x_post;      // ambient temperature
  vector[N_post] y_post;      // outcome energy vector
}
parameters {
  // This block declares the parameters of the model.
  real alpha;      // baseline consumption
  real beta_h;     // slope for heating
  real tau_h;      // threshold temperature for heating
  real beta_c;     // slope for cooling
  real tau_c;      // threshold temperature for cooling
  real<lower=0> sigma;  // error scale
  real theta;     // moving average
}
transformed parameters {
  vector[N_pre] f_pre;
  vector[N_pre] epsilon;
  vector[N_post] f_post;
  for (n in 1:N_pre) {
    f_pre[n] = alpha + beta_h * fmax(tau_h-x_pre[n], 0) + beta_c * fmax(x_pre[n]-tau_c, 0);
  }
  epsilon[1] = y_pre[1] - f_pre[1];
  for (n in 2:N_pre) {
    epsilon[n] = y_pre[n] - f_pre[n] - theta*epsilon[n-1];
  }
  for (n in 1:N_post) {
    f_post[n] = alpha + beta_h * fmax(tau_h-x_post[n], 0) + beta_c * fmax(x_post[n]-tau_c, 0);
  }
}
model {
  // Assigning prior distributions on some parameters
  alpha ~ normal(800, 50);
  tau_h ~ normal(6, 2);
  tau_c ~ normal(16, 2);
  beta_h ~ normal(30, 5);
  beta_c ~ normal(30, 5);
  theta ~ normal(0.7, 0.2);
  sigma ~ normal(80, 10);
  // Observational model
  for (n in 2:N_pre) {
    y_pre[n] ~ normal(f_pre[n] + theta*epsilon[n-1], sigma);
  }
}
generated quantities {
  // This block is for posterior predictions. It is not part of model training
  vector[N_pre] y_pre_pred;
  vector[N_post] y_post_pred;
  real savings = 0;
  
  y_pre_pred[1] = normal_rng(f_pre[1], sigma);
  for (n in 2:N_pre) {
    y_pre_pred[n] = normal_rng(f_pre[n] + theta*epsilon[n-1], sigma);
  }
  
  // Forecasts with increased uncertainty due to the autocorrelation
  y_post_pred[1] = normal_rng(f_post[1], sigma);
  for (n in 2:N_post) {
    y_post_pred[n] = normal_rng(f_post[n], sigma * sqrt(1+theta^2));
    savings += y_post_pred[n] - y_post[n];
  }
}
