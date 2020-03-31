data {
  int<lower=0> N;   // number of trials
  int<lower=0> K;   // number of time points
  matrix[K, N] X;   // predictor matrix (differences)
  int<lower=0,upper=1> Y[N];      // outcome vector
}
parameters {
  real alpha;                 // location (bias)
  row_vector[K] beta_raw;     // coefficients for predictors
  real<lower=0> sigma;        // scale
  real<lower=0> tau;          // smoothness
}
transformed parameters{
  row_vector[N] mu;
  row_vector[K] beta;         // coefficient vector - random walk prior
  beta[K] = beta_raw[K];
  for (i in 2:K){
    beta[K+1-i] = beta[K+2-i] + (beta_raw[i]*tau); 
  }
  mu = (alpha + beta*X) ./(sqrt(2)*sigma);
}
model {
  // priors
  sigma ~ cauchy(0, 1);
  alpha ~ normal(0, 1);
  beta_raw ~ normal(0, 0.1);
  tau ~ normal(0, 0.1);
  
  // likelihood
  for (n in 1:N){
    Y[n] ~ bernoulli( Phi(mu[n]) );  
  }
}
