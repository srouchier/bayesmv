data {
  // This block declares all data which will be passed to the Stan model.
  // Pre period
  int<lower=0> N_pre;       // number of data items
  vector[N_pre] x_pre;      // outdoor temperature
  vector[N_pre] y_pre;      // outcome energy vector
  // Post period
  int<lower=0> N_post;        // number of data items
  vector[N_post] x_post;      // outdoor temperature
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
}
model {
  // Assigning prior distributions on some parameters
  alpha ~ normal(800, 100);
  tau_h ~ normal(8, 5);
  tau_c ~ normal(18, 5);
  beta_h ~ normal(40, 15);
  beta_c ~ normal(40, 15);
  // Observational model
  for (n in 1:N_pre) {
    y_pre[n] ~ normal(alpha + beta_h * fmax(tau_h-x_pre[n], 0) + beta_c * fmax(x_pre[n]-tau_c, 0), sigma);
  }
}
generated quantities {
  // This block is for posterior predictions. It is not part of model training
  vector[N_pre] y_pre_pred;
  vector[N_post] y_post_pred;
  real savings = 0;
  
  for (n in 1:N_pre) {
    y_pre_pred[n] = normal_rng(alpha + beta_h * fmax(tau_h-x_pre[n], 0) + beta_c * fmax(x_pre[n]-tau_c, 0), sigma);
  }
  
  for (n in 1:N_post) {
    y_post_pred[n] = normal_rng(alpha + beta_h * fmax(tau_h-x_post[n], 0) + beta_c * fmax(x_post[n]-tau_c, 0), sigma);
    savings += y_post_pred[n] - y_post[n];
  }
}
