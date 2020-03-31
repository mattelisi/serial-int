# --------------------------------------------------------------- #
# This script prepare the data for the modelling in Stan
#
# In the main dataset (made with 'build_dataset.R') there are the onset times
# of each stimulus sample relative to the saccade onset. The code below expand the 
# the stimulus and calculate the average position/luminance of the stimulus along
# 100 equally spaced temporal intervals. Within any given time bin the value of 
# the stimulus is obtained by averaging the stimulus value, weighting by their
# relative prevalence in the interval.
#
# Matteo Lisi 2019
# --------------------------------------------------------------- #


rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")
system("mkdir stan_data_xp1") # make folder

# load dataset
d_all <- read.table("./data/exp1_data",sep="\t",header=T)

for(id_i in 1:length(unique(d_all$id))){
  
  id_label <- unique(d_all$id)[id_i]
  d <- d_all[d_all$id==id_label,]
  
  trial_label <- paste(d$vpcode, d$vp,d$session,d$block,d$trial3,sep="_")
  d$trial <- as.numeric(factor(trial_label, labels=1:length(unique(trial_label))))
  
  # build data matrix
  range(d$tSteps)
  n_time_bin <- 100
  bin_width <- 900/n_time_bin
  
  # saccade data
  X_h <- matrix(NA,nrow=100, ncol=max(d$trial))
  X_v <- matrix(NA,nrow=100, ncol=max(d$trial))
  S_h <- rep(NA, max(d$trial))
  S_v <- rep(NA, max(d$trial))
  
  # perception data
  PosDiff <- matrix(NA,nrow=100, ncol=max(d$trial))
  choice <- rep(NA, max(d$trial))
  
  # populate
  for(i in unique(d$trial)){
    D_di <- d$tSteps[d$trial==i]
    x_di <- d$xs[d$trial==i] + d$tarX[d$trial==i]
    y_di <- d$ys[d$trial==i]
    S_h[i] <- unique(d$sacXresp[d$trial==i]) # maybe add sanity check that this is a scalar?
    S_v[i] <- unique(d$sacYresp[d$trial==i])
    
    di_xR <- ifelse(d$side[d$trial==i]==1, d$xNear[d$trial==i], d$xFar[d$trial==i])
    di_xL <- ifelse(d$side[d$trial==i]==-1, d$xNear[d$trial==i], d$xFar[d$trial==i])
    di_yR <- ifelse(d$side[d$trial==i]==1, d$yNear[d$trial==i],d$yFar[d$trial==i])
    di_yL <- ifelse(d$side[d$trial==i]==-1, d$yNear[d$trial==i], d$yFar[d$trial==i])
    di_xyR <- sqrt(di_yR^2 + di_xR^2)
    di_xyL <- sqrt(di_yL^2 + di_xL^2)
    di_posDiff <- di_xyR - di_xyL
    choice[i] <- 1-unique(d$resp[d$trial==i])
    
    # expand
    D_i <- seq(-900,0)
    x_i <- rep(0, length(D_i))
    y_i <-  rep(0, length(D_i))
    PD_i <- rep(0, length(D_i))
    current <- 1
    for(t_i in 1:length(D_i)){
      if(D_i[t_i] >= D_di[current]){
        current <- current + 1
      }
      if(current > 1){
        if(D_i[t_i]>=D_di[current-1]){
          x_i[t_i] <- x_di[current-1]
          y_i[t_i] <- y_di[current-1]
          PD_i[t_i] <- di_posDiff[current-1]
        }
      }
    }
    
    # # sanity check
    # plot(D_i, x_i, type="l") # sanity check
    # points(D_di, x_di)
    
    #this is discrete bin
    for(t_i in 1:n_time_bin){
      index_t <- (1+(t_i-1)*bin_width):(t_i*bin_width)
      X_h[t_i,i] <- mean(x_i[index_t],na.rm=T)
      X_v[t_i,i] <- mean(y_i[index_t],na.rm=T)
      PosDiff[t_i,i] <- mean(PD_i[index_t],na.rm=T)
    }
    
  }
  
  # make data list for fitting with Stan
  # saccade
  d_stan <- list(Xh = X_h, Xv = X_v, Yh= S_h, Yv=S_v, K=n_time_bin, N=max(d$trial))
  str(d_stan)
  # perception
  d_prc <- list(X = PosDiff, Y = choice, K=n_time_bin, N=max(d$trial))
  str(d_prc)
  
  saveRDS(d_stan, paste("./stan_data_xp1/",id_label,"_saccade_xp1.RDS",sep=""))
  saveRDS(d_prc, paste("./stan_data_xp1/",id_label,"_perception_xp1.RDS",sep=""))
}


# --------------------------------------------------------------- #
# same but for exp 2
rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")
system("mkdir stan_data_xp2")

# load dataset
d_all <- read.table("./data/exp2_data",sep="\t",header=T)

for(id_i in 1:length(unique(d_all$id))){
  
  id_label <- unique(d_all$id)[id_i]
  d <- d_all[d_all$id==id_label,]
  
  trial_label <- paste(d$vpcode, d$vp,d$session,d$block,d$trial3,sep="_")
  d$trial <- as.numeric(factor(trial_label, labels=1:length(unique(trial_label))))
  
  # build data matrix
  range(d$tSteps)
  n_time_bin <- 100
  bin_width <- 900/n_time_bin
  
  # saccade data
  X_h <- matrix(NA,nrow=100, ncol=max(d$trial))
  X_v <- matrix(NA,nrow=100, ncol=max(d$trial))
  S_h <- rep(NA, max(d$trial))
  S_v <- rep(NA, max(d$trial))
  
  # perception data
  LumDiff <- matrix(NA,nrow=100, ncol=max(d$trial))
  choice <- rep(NA, max(d$trial))
  
  # populate
  for(i in unique(d$trial)){
    D_di <- d$tSteps[d$trial==i]
    x_di <- d$xs[d$trial==i] + d$tarX[d$trial==i]
    y_di <- d$ys[d$trial==i]
    S_h[i] <- unique(d$sacXresp[d$trial==i]) # maybe add sanity check that this is a scalar?
    S_v[i] <- unique(d$sacYresp[d$trial==i])
    
    lum_di <- d$LumDiff[d$trial==i]
    choice[i] <- 1-unique(d$resp[d$trial==i])
    
    # expand
    D_i <- seq(-900,0)
    x_i <- rep(0, length(D_i))
    y_i <-  rep(0, length(D_i))
    lum_i <- rep(0, length(D_i))
    current <- 1
    for(t_i in 1:length(D_i)){
      if(D_i[t_i] >= D_di[current]){
        current <- current + 1
      }
      if(current > 1){
        if(D_i[t_i]>=D_di[current-1]){
          x_i[t_i] <- x_di[current-1]
          y_i[t_i] <- y_di[current-1]
          lum_i[t_i] <- lum_di[current-1]
        }
      }
    }
    
    # # sanity check
    # plot(D_i, x_i, type="l") # sanity check
    # points(D_di, x_di)
    
    #this is discrete bin
    for(t_i in 1:n_time_bin){
      index_t <- (1+(t_i-1)*bin_width):(t_i*bin_width)
      X_h[t_i,i] <- mean(x_i[index_t],na.rm=T)
      X_v[t_i,i] <- mean(y_i[index_t],na.rm=T)
      LumDiff[t_i,i] <- mean(lum_i[index_t],na.rm=T)
    }
    
  }
  
  # make data list for fitting with Stan
  # saccade
  d_stan <- list(Xh = X_h, Xv = X_v, Yh= S_h, Yv=S_v, K=n_time_bin, N=max(d$trial))
  str(d_stan)
  # perception
  d_prc <- list(X = LumDiff, Y = choice, K=n_time_bin, N=max(d$trial))
  str(d_prc)
  
  saveRDS(d_stan, paste("./stan_data_xp2/",id_label,"_saccade_xp2.RDS",sep=""))
  saveRDS(d_prc, paste("./stan_data_xp2/",id_label,"_perception_xp2.RDS",sep=""))
}
