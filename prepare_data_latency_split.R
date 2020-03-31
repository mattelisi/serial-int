# --------------------------------------------------------------- #
# This script prepare the data for the modelling in Stan
#
# Latency split analysis
#

# --------------------------------------------------------------- #
rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")
rm(list=ls())

# --------------------------------------------------------------- #
# dataset exp 1
d <- read.table("./data/exp1_data",sep="\t",header=T)
d$id <- paste(d$id,"_xp1",sep="")

dxy <- data.frame(D = d$tSteps,  sacXresp=d$sacXresp, sacYresp= d$sacYresp, xs = d$xs, ys = d$ys, meanXpos = abs(d$tarX), id=d$id, rt=d$sacOnset, vpcode=d$vpcode, block=d$block, trial3=d$trial3,session=d$session)

d$xR <- ifelse(d$side==1, d$xNear, d$xFar)
d$xL <- ifelse(d$side==-1, d$xNear, d$xFar)
d$yR <- ifelse(d$side==1, d$yNear,d$yFar)
d$yL <- ifelse(d$side==-1, d$yNear, d$yFar)
d$xyR <- sqrt(d$yR^2 + d$xR^2)
d$xyL <- sqrt(d$yL^2 + d$xL^2)

d$xRmean <- ifelse(d$side==1, d$meanXnear, d$meanXfar)
d$xLmean <- ifelse(d$side==-1, d$meanXnear, d$meanXfar)

dpe <- data.frame(D = d$tSteps, resp=-d$resp+1, xDiff =  d$xyR-d$xyL, xDiffMean = d$xRmean-abs(d$xLmean), vpcode=d$vpcode, block=d$block, trial3=d$trial3, session=d$session, id=d$id, rt=d$sacOnset)

dxy_1 <- dxy
dpe_1 <- dpe
rm(dxy, dpe)

# --------------------------------------------------------------- #
# dataset exp 2
d <- read.table("./data/exp2_data",sep="\t",header=T)
d$id <- paste(d$id,"_xp2",sep="")

dxy <- data.frame(D = d$tSteps, sacXresp=d$sacXresp, sacYresp= d$sacYresp, xs = d$xs, ys = d$ys, meanXpos = abs(d$tarX), id=d$id, rt=d$sacOnset, vpcode=d$vpcode, block=d$block, trial3=d$trial3,session=d$session)

dpe <- data.frame(D = d$tSteps, resp=-d$resp+1, xDiff = d$LumDiff, xDiffMean = d$meanLumDiff, vpcode=d$vpcode, block=d$block, trial3=d$trial3, session=d$session, id=d$id, rt=d$sacOnset)

dxy_2 <- dxy
dpe_2 <- dpe
rm(dxy, dpe)

# --------------------------------------------------------------- #
# equalize values of perceptual task stimuli
str(dpe_1)
str(dpe_2)
dpe_1$xDiff <- scale(dpe_1$xDiff, center=F) * 2.5 # just to have perceptual stimuli in the same scale 
dpe_2$xDiff <- scale(dpe_2$xDiff, center=F) * 2.5
dpe_1$xDiffMean <- scale(dpe_1$xDiffMean, center=F) * 2.5
dpe_2$xDiffMean <- scale(dpe_2$xDiffMean, center=F) * 2.5

# --------------------------------------------------------------- #
# merge database
dpe <- rbind(dpe_1, dpe_2)
dxy <- rbind(dxy_1, dxy_2)
rm(dpe_1,dpe_2,dxy_1,dxy_2)

# --------------------------------------------------------------- #
# compute latency bins
dpe$bin <- NA
dxy$bin <- NA
for(i in unique(dpe$id)){
  # quartiles
  dpe$bin[dpe$id==i] <- cut(dpe$rt[dpe$id==i],breaks=quantile(dpe$rt[dpe$id==i]),labels=1:4,include.lowest=T)
  dxy$bin[dxy$id==i] <- cut(dxy$rt[dxy$id==i],breaks=quantile(dxy$rt[dxy$id==i]),labels=1:4,include.lowest=T)
}

any(is.na(dpe$bin))
any(is.na(dxy$bin))

## 
saveRDS(dxy, paste("./data/rt_split/raw_data_saccade_RTbin.RDS"))
saveRDS(dpe, paste("./data/rt_split/raw_data_perception_RTbin.RDS"))


unique(dxy$bin)

# --------------------------------------------------------------- #
# make stan data for each bin

# trial_label <- paste(dpe$vpcode, dpe$vp,dpe$session,dpe$block,dpe$trial3,sep="_")
# dpe$trial <- as.numeric(factor(trial_label, labels=1:length(unique(trial_label))))
# count_uniques <- function(x){ length(unique(x))}
# trial_count_table <- tapply(dpe$trial, list(dpe$bin, dpe$id), count_uniques)

for(bin_i in 1:4){
  
  for(id_i in 1:length(unique(dpe$id))){
    
    id_label <- unique(dpe$id)[id_i]
    d_sac_i <- dxy[dxy$id==id_label & dxy$bin==bin_i,]
    d_per_i <- dpe[dpe$id==id_label & dpe$bin==bin_i,]
    
    trial_label <- paste(d_sac_i$vpcode, d_sac_i$vp,d_sac_i$session,d_sac_i$block,d_sac_i$trial3,sep="_")
    d_sac_i$trial <- as.numeric(factor(trial_label, labels=1:length(unique(trial_label))))
    d_per_i$trial <- as.numeric(factor(trial_label, labels=1:length(unique(trial_label))))
    
    # build data matrix
    # range(d_sac_i$D)
    n_time_bin <- 100
    bin_width <- 900/n_time_bin
    
    # saccade data
    X_h <- matrix(NA,nrow=100, ncol=max(d_sac_i$trial))
    X_v <- matrix(NA,nrow=100, ncol=max(d_sac_i$trial))
    S_h <- rep(NA, max(d_sac_i$trial))
    S_v <- rep(NA, max(d_sac_i$trial))
    
    # perception data
    xDiff <- matrix(NA,nrow=100, ncol=max(d_sac_i$trial))
    choice <- rep(NA, max(d_sac_i$trial))
    
    # populate
    for(i in sort(unique(d_sac_i$trial))){
      D_di <- d_sac_i$D[d_sac_i$trial==i]
      x_di <- d_sac_i$xs[d_sac_i$trial==i] + d_sac_i$meanXpos[d_sac_i$trial==i]
      y_di <- d_sac_i$ys[d_sac_i$trial==i]
      S_h[i] <- unique(d_sac_i$sacXresp[d_sac_i$trial==i]) # maybe add sanity check that this is a scalar?
      S_v[i] <- unique(d_sac_i$sacYresp[d_sac_i$trial==i])
      
      xDiff_di <- d_per_i$xDiff[d_per_i$trial==i]
      # choice[i, id_i] <- 1-unique(d_per_i$resp[d_per_i$trial==i])
      choice[i] <- unique(d_per_i$resp[d_per_i$trial==i]) # already inverted above
      
      # expand
      D_i <- seq(-900,0)
      x_i <- rep(0, length(D_i))
      y_i <-  rep(0, length(D_i))
      xDiff_i <- rep(0, length(D_i))
      current <- 1
      for(t_i in 1:length(D_i)){
        if(D_i[t_i] >= D_di[current]){
          current <- current + 1
        }
        if(current > 1){
          if(D_i[t_i]>=D_di[current-1]){
            x_i[t_i] <- x_di[current-1]
            y_i[t_i] <- y_di[current-1]
            xDiff_i[t_i] <- xDiff_di[current-1]
          }
        }
      }
      
      # # sanity check
      # plot(D_i, x_i, type="l") # sanity check
      # points(D_di, x_di)
      
      #this is discrete bin
      for(t_i in 1:n_time_bin){
        index_t <- (1+(t_i-1)*bin_width):(t_i*bin_width)
        X_h[ t_i, i] <- mean(x_i[index_t],na.rm=T)
        X_v[t_i, i] <- mean(y_i[index_t],na.rm=T)
        xDiff[ t_i, i] <- mean(xDiff_i[index_t],na.rm=T)
      }
    }
    
    # make data list for fitting with Stan
    # saccade
    d_stan <- list(Xh = X_h, Xv = X_v, Yh= S_h, Yv=S_v, K=n_time_bin, N=max(d_sac_i$trial))
    str(d_stan)
    # perception
    d_prc <- list(X = xDiff, Y = choice, K=n_time_bin, N=max(d_sac_i$trial))
    str(d_prc)
    
    saveRDS(d_stan, paste("./stan_data_rtsplit/",id_label,"_bin",bin_i,"saccade.RDS",sep=""))
    saveRDS(d_prc, paste("./stan_data_rtsplit/",id_label,"_bin",bin_i,"perception.RDS",sep=""))
    
  }
}
