# --------------------------------------------------------------- #
# run all analyses for Exp. 2

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")
system("mkdir ./results/xp2/")

# load dataset
d_all <- read.table("./data/exp2_data",sep="\t",header=T)
ID_s <- unique(d_all$id)
rm(d_all)

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

# run stan models, looping over subjects
for(id_i in ID_s){
  d_stan <- readRDS(paste("./stan_data_xp2/",id_i,"_saccade_xp2.RDS",sep=""))
  m0 <- stan(file = "./stan-code/saccade_model.stan", data = d_stan, iter = 3000, chains = 4, control = list(max_treedepth = 12))
  saveRDS(m0,paste("./results/xp2/",id_i,"_saccade_xp2_fit.RDS",sep=""))
  rm(m0, d_stan)
  
  d_prc <- readRDS(paste("./stan_data_xp2/",id_i,"_perception_xp2.RDS",sep=""))
  m1 <- stan(file = "./stan-code/perception_model.stan", data = d_prc, iter = 3000, chains = 4, init=init_per )
  saveRDS(m1,paste("./results/xp2/",id_i,"_perception_xp2_fit.RDS",sep=""))
  rm(m1, d_prc)
}
