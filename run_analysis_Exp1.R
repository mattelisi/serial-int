# --------------------------------------------------------------- #
# run all analyses for exp. 1

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/") # you may have to change the folder
system("mkdir ./results/xp1/")

# load dataset (for ID list)
d_all <- read.table("./data/exp1_data",sep="\t",header=T)
ID_s <- unique(d_all$id)
rm(d_all)

# run model
library(rstan)
options(mc.cores = parallel::detectCores()) # use all cores

# function to set initial values in perception model
init_per <- function() {
  list(alpha = rnorm(1,0,0.2), 
       sigma = exp(rnorm(1,0.5,0.2)), 
       beta_raw = rnorm(100,0,0.01), 
       tau = exp(rnorm(1,-1,0.05)),
       beta=rnorm(100,0,0.01))
}

# run stan models, looping over subjects
for(id_i in ID_s){
  d_stan <- readRDS(paste("./stan_data_xp1/",id_i,"_saccade_xp1.RDS",sep=""))
  m0 <- stan(file = "./stan-code/saccade_model.stan", data = d_stan, iter = 3000, chains = 4, control = list(max_treedepth = 12))
  saveRDS(m0,paste("./results/xp1/",id_i,"_saccade_xp1_fit.RDS",sep=""))
  rm(m0, d_stan)
  
  d_prc <- readRDS(paste("./stan_data_xp1/",id_i,"_perception_xp1.RDS",sep=""))
  m1 <- stan(file = "./stan-code/perception_model.stan", data = d_prc, iter = 3000, chains = 4, init=init_per )
  saveRDS(m1,paste("./results/xp1/",id_i,"_perception_xp1_fit.RDS",sep=""))
  rm(m1, d_prc)
}
