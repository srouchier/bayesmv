data {
// Baseline data: pre ECM
  int<lower=0> N;     // number of data items
  int<lower=0> K;     // number of predictors
  matrix[N, K] x;     // predictor matrix
  vector[N] y;        // outcome vector
// Reporting period data: post ECM
  int<lower=0> N_post;  // number of data items
  matrix[N, K] x_post;     // predictor matrix
  vector[N] y_post;        // outcome vector
}

parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
 }
 
model {
  y ~ normal(x * beta + alpha, sigma);  // likelihood
}

generated quantities {
  vector[N_post] prediction;
  real savings = 0;
  for (n in 1:N_post) {
    prediction[n] = normal_rng(x_post[n] * beta + alpha, sigma);
    savings += prediction[n] - y_post[n];
  }
}