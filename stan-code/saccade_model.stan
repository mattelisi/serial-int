data {
  int<lower=0> N;   // number of trials
  int<lower=0> K;   // number of time points
  matrix[K, N] Xh;  // predictor matrix (horiz.)
  matrix[K, N] Xv;  // predictor matrix (vertical.)
  vector[N] Yh;     // outcome vector (horiz.)
  vector[N] Yv;     // outcome vector (vert.)
}
parameters {
  vector[2] alpha;            // intercept (bias)
  real vScale;                // scaling for vertical coordinates
  row_vector[K] beta_raw;     // coefficients for predictors
  vector<lower=0>[2] sigma;   // error scale
  real<lower=0> tau;          // smoothness
}
transformed parameters{
  row_vector[K] beta;         // coefficient vector - random walk prior
  beta[K] = beta_raw[K];
  for (i in 2:K){
    beta[K+1-i] = beta[K+2-i] + (beta_raw[i]*tau); 
  }
}
model {
  // priors
  sigma ~ cauchy(0, 1);
  alpha ~ normal(0, 1);
  vScale ~ normal(1, 1);
  beta_raw ~ normal(0, 0.1);
  tau ~ normal(0, 0.1);
  
  // likelihood
  Yh ~ normal(alpha[1] + beta*Xh, sigma[1]);  
  Yv ~ normal(alpha[2] + (vScale*beta)*Xv, sigma[2]);  
}
