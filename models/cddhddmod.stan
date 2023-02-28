data {
  int<lower=0> N;     // number of data items
  int<lower=0> K;     // number of predictors
  matrix[N, K] x;     // predictor matrix
  vector[N] y;        // outcome vector
  real cddpost;       // post-retrofit cdd
  real hddpost;       // post-retrofit hdd
  real enrpost;       // post-retrofit annual energy use
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
  real n12;             // annual error
 }
model {
  y ~ normal(x * beta + alpha, sigma);  // likelihood
  n12 ~ normal(0,sqrt(12.0)*sigma);
}
generated quantities {
  real svgs;
  svgs = beta[1] * cddpost + beta[2]*hddpost + 12.0*alpha + n12 - enrpost;
}
