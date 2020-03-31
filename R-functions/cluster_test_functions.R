# Set of custom functions to perform 1D cluster test.
# The functions calculate cluster probability for both Z and t fields
# Matteo Lisi, 2019

## expected Euler characteristics
expected_euler_1D <- function(th){
  (sqrt(4 * log(2))/(2*pi)) * exp(-(th^2)/2)
}
# for a t field (v are the degree of freedom)
expected_euler_1D_t <- function(th, v){
  (sqrt(4 * log(2))/(2*pi)) * (1 + (th^2)/v)^(-(v-1)/2)
}

## probability of cluster (for z and t fields)
# th: threshold used
# k:  size of cluster
# S:  total size (length, since 1D) of search space
# R:  number of resels
cluster_p <- function(th, k, S, R){
  E_m <- R*expected_euler_1D(th) # expectred Euler C. for threshold th
  E_N <- S * pnorm(-th)
  beta <- (gamma(3/2) * (E_m/E_N))^2
  p_cluster <- 1 - exp(-E_m * exp(-beta * k^2))
  return(p_cluster)
}
# for a t field (v are the degree of freedom)
cluster_p_t <- function(th, k, S, R, v){
  E_m <- R*expected_euler_1D_t(th, v)
  E_N <- S * pt(-th, df=v)
  beta <- (gamma(3/2) * (E_m/E_N))^2
  p_cluster <- 1 - exp(-E_m * exp(-beta * k^2))
  return(p_cluster)
}
