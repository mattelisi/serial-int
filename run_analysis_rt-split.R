# run latency-split analysis

# individual fits
rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")
system("mkdir ./results/rt-split/") # make folder

# run model
library(rstan)
options(mc.cores = parallel::detectCores())

# function to set initial values in perception model
init_per <- function() {
  list(alpha = rnorm(1,0,0.2), 
       sigma = exp(rnorm(1,0.5,0.2)), 
       beta_raw = rnorm(100,0,0.001), 
       tau = exp(rnorm(1,-1,0.05)),
       beta=rnorm(100,0,0.001))
}

# get list of id_s
dxy <- readRDS("./data/rt_split/raw_data_saccade_RTbin.RDS")
ID_s <- as.character(unique(dxy$id))
rm(dxy)

# run stan models, looping over subjects
for(bin_i in 1:4){
  for(id_i in ID_s){
    
    d_stan <- readRDS(paste("./stan_data_rtsplit/",id_i,"_bin",bin_i,"saccade.RDS",sep=""))
    m0 <- stan(file = "./stan-code/saccade_model.stan", data = d_stan, iter = 3000, chains = 4, control = list(max_treedepth = 10))
    saveRDS(m0,paste("./results/rt-split/",id_i,"_bin",bin_i,"_fit.RDS",sep=""))
    rm(m0, d_stan)
    
    d_prc <- readRDS(paste("./stan_data_rtsplit/",id_i,"_bin",bin_i,"perception.RDS",sep=""))
    m1 <- stan(file = "./stan-code/perception_model.stan", data = d_prc, iter = 3000, chains = 4, init=init_per )
    saveRDS(m1,paste("./results/rt-split/",id_i,"_bin",bin_i,"_fitPerception.RDS",sep=""))
    rm(m1, d_prc)
    
  }
}
