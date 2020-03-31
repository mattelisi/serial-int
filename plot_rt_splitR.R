# --------------------------------------------------------------- #
# plot group level functions
# Bayesian (onestep), RT split analysis

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")

# get list of id_s
dxy <- readRDS("./data/rt_split/raw_data_saccade_RTbin.RDS")
ID_s <- unique(dxy$id)
rm(dxy)

library(rstan)
library(rethinking)
library(ggplot2)

# this is only for loading a custom ggplot theme...
source("/home/matteo/sync/miscR/miscFunctions.R") 

# here is the theme I used:
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

# load individual functions
d <- {}
for(bin_i in 1:4){
  for(id_i in ID_s){
    
    file_saccade <- paste("./results/rt-split/",id_i,"_bin",bin_i,"_fit.RDS",sep="")
    file_perception <- paste("./results/rt-split/",id_i,"_bin",bin_i,"_fitPerception.RDS",sep="")

    m0 <- readRDS(file_saccade)
    m1 <- readRDS(file_perception)
    
    beta_samples <- extract(m0, pars=c("beta"))$beta
    plo <- data.frame(beta=apply(beta_samples, 2, mean))
    plo$SE <- data.frame(beta=apply(beta_samples, 2, sd))$beta
    plo$lb <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[1])
    plo$ub <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[2])
    plo$lb_se <- plo$beta - apply(beta_samples, 2, sd)
    plo$ub_se <- plo$beta + apply(beta_samples, 2, sd)
    plo$D <- seq(-900,-9,length.out=100)+9/2 # bin centers
    plo$signif <- ifelse(plo$lb>0 | plo$ub<0,1,0)
    plo$signif <- ifelse(plo$signif, plo$signif * (-0.04), NA)
    plo$type <- "saccade"
    
    beta2_samples <- extract(m1, pars=c("beta"))$beta 
    plo2 <- data.frame(beta=apply(beta2_samples, 2, mean))
    plo2$SE <- data.frame(beta=apply(beta2_samples, 2, sd))$beta
    plo2$lb <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[1])
    plo2$ub <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[2])
    plo2$lb_se <- plo2$beta - apply(beta2_samples, 2, sd)
    plo2$ub_se <- plo2$beta + apply(beta2_samples, 2, sd)
    plo2$D <- seq(-900,-9,length.out=100)+9/2
    plo2$signif <- ifelse(plo2$lb>0 | plo2$ub<0,1,0)
    plo2$signif <- ifelse(plo2$signif, plo2$signif * (-0.04), NA)
    plo2$type <- "perception"
    
    plo <- rbind(plo,plo2)
    
    
    plo$id <- id_i
    plo$bin <- bin_i
    d <- rbind(d, plo)
    
  }
}

d$signi_plot <-ifelse(!is.na(d$signif)& d$type=="saccade",-0.065, d$signif)

# --------------------------------------------------------------- #
# make group-level plot

colnames(d)[1] <- "w"
N <- length(unique(d$id))
group_dat <- d

# average and compute standard error across participants
dag <- aggregate(w~D+type+bin,group_dat, mean)
dag$se <- aggregate(w~D+type+bin,group_dat, se)$w
str(dag)

# set a probability threshold
alpha <- 0.05

# thresholding
group_dat$signif <- NA
group_dat$signi_plot <- NA
group_dat$signi <- ifelse(group_dat$lb>0 | group_dat$ub<0,1,NA)

#cluster test - custom functions
source("./R-functions/cluster_test_functions.R")

## set parameters for cluster test
n_resels <- 4*900 / (1/15 * 1000)
th <- qnorm(1-alpha/2)
p_threhold_cluster <- 0.01

inblob <- 0; start_blob <- NA; end_blob <- NA; blob_p<-NA
start_blob_D <- NA; end_blob_D <- NA;
count <- 0; n_blob <- 0; d_blob_out <- {}
for(i in 2:length(group_dat$signi)){
  if(!is.na(group_dat$signi[i]) & is.na(group_dat$signi[i-1])){
    inblob <- 1
    start_blob <- i
    start_blob_D <- group_dat$D[i]
  }
  if(!is.na(group_dat$signi[i-1]) & is.na(group_dat$signi[i]) & inblob==1){
    inblob <- 0
    end_blob <- i
    end_blob_D <- group_dat$D[i]
    blob_length <- end_blob - start_blob
    n_blob <-n_blob +1
    
    blob_p <- cluster_p(th, k=blob_length, S=400, R=n_resels)
    
    d_line <- data.frame(n=n_blob, s=blob_length, type=group_dat$type[i-1],p = blob_p, blob_start=start_blob, blob_end = end_blob, blob_start_D= start_blob_D, blob_end_D = end_blob_D, id=group_dat$id[i], bin=group_dat$bin[i])
    d_blob_out <- rbind(d_blob_out, d_line)
    
    if(blob_p>=p_threhold_cluster){
      group_dat$signi[start_blob:end_blob]<-NA
      count = count+1
      cat("- removed blob ",count," of size ",blob_length,"\n")
    }
  }
}
d_blob_out$excl <- ifelse(d_blob_out$p >= p_threhold_cluster, 1,0)
d_blob_out <- d_blob_out[d_blob_out$excl==0,]

# change based on ID the 'height' of signi_plot
# used to add "significance" column to the aggregate database
# based on average onset-offset of the integration window 
group_dat$id_n <- as.numeric(as.factor(group_dat$id))
group_dat$signi_plot <- group_dat$signi* (-0.0135) - 0.004*group_dat$id_n
group_dat$se <- NA

# calculate the mean start and end points of the intergration windows
# for the saccade there is only 1 blob for each participants here so it's OK
# just first remove early spurious saccade clusters (e.g. 800 msec before), if present
d_blob_out <- d_blob_out[-which(d_blob_out$blob_end_D< (-800) & d_blob_out$type=="saccade"),]
saccade_avg <- aggregate(cbind(blob_start_D,blob_end_D)~bin, d_blob_out[d_blob_out$type=="saccade",], mean)
saccade_avg$start_se <- aggregate(blob_start_D~bin, d_blob_out[d_blob_out$type=="saccade",], se)$blob_start_D
saccade_avg$end_se <- aggregate(blob_end_D~bin, d_blob_out[d_blob_out$type=="saccade",], se)$blob_end_D
saccade_avg[,2:5] <- saccade_avg[,2:5]/1000

# for the perception I take and align only the bigger blob
blob_per <- d_blob_out[d_blob_out$type=="perception",]

# remove the smaller early blob for 1 subject (only for the averaging)
Ntab <- tapply(blob_per$id, list(blob_per$id, blob_per$bin), length)
id_extra <- names(which(apply(Ntab,1,sum,na.rm=T)>4))
bin_extra <- which(Ntab[which(apply(Ntab,1,sum,na.rm=T)>4),]>1)
blob_sizes <- blob_per$s[blob_per$id==id_extra & blob_per$bin==bin_extra]
blob_per <- blob_per[-which(blob_per$id==id_extra & blob_per$bin==bin_extra & blob_per$s==min(blob_sizes)),]

perception_avg <- aggregate(cbind(blob_start_D,blob_end_D)~bin, blob_per, mean)
perception_avg$start_se <- aggregate(blob_start_D~bin, blob_per, se)$blob_start_D
perception_avg$end_se <- aggregate(blob_end_D~bin, blob_per, se)$blob_end_D
perception_avg[,2:5] <- perception_avg[,2:5]/1000

# group y coordinate for segments
y_group <- (-0.02) - 0.004*(max(group_dat$id_n)+1)

# prepare the dataframe for plotting
group_dat$id <- factor(group_dat$id)
group_dat$type <- factor(group_dat$type)
saccade_avg$y_fixed <- y_group
perception_avg$y_fixed <- y_group

# finally make plot
require("gridExtra") 

#group_plot <- 
# dag,aes(x=D/1000,y=w,ymin=w-se,ymax=w+se,color=type, fill=type)
kernel_plot <- (ggplot() +geom_line(data=group_dat,aes(x=D/1000,y=w, group=id:type, color=type),size=0.4,alpha=0.1)
                +geom_ribbon(data=dag, alpha=0.25,color=NA, aes(x=D/1000,ymin=w-se,ymax=w+se, fill=type))
                +geom_line(data=dag,aes(x=D/1000,y=w,color=type),size=1.2)
                +geom_hline(yintercept=0,lty=2)
                +coord_cartesian(ylim=c(-0.062,0.12),xlim=c(-950,-38)/1000)
                +geom_vline(xintercept=0,lty=1)
                +scale_x_continuous(breaks=seq(-1000,0,200)/1000)
                +scale_y_continuous(breaks=seq(0,1,0.05))
                +nice_theme+scale_color_manual(values=c("black","blue"),name="",guide=F)
                +scale_fill_manual(values=c("black","blue"),name="",guide=F)
                +facet_grid(.~bin)
                +geom_line(data=group_dat, aes(x=D/1000,y=signi_plot,group=id:type, color=type),size=0.9)
                +labs(x="time relative to saccade onset [sec]",y="weight")
                +geom_segment(data= saccade_avg, aes(x=blob_start_D,y=y_fixed, xend=blob_end_D,yend=y_fixed),size=5, color="blue", alpha=0.3)
                +geom_segment(data= saccade_avg, aes(x=blob_start_D-start_se, y=y_fixed, xend=blob_start_D+start_se, yend=y_fixed),color="blue")
                +geom_segment(data= saccade_avg, aes(x=blob_end_D-end_se, y=y_fixed, xend=blob_end_D+end_se, yend=y_fixed),color="blue")
                +geom_segment(data= perception_avg, aes(x=blob_start_D,y=y_fixed, xend=blob_end_D,yend=y_fixed),size=5, color="black", alpha=0.3)
                +geom_segment(data= perception_avg, aes(x=blob_start_D-start_se, y=y_fixed, xend=blob_start_D+start_se, yend=y_fixed),color="black")
                +geom_segment(data= perception_avg, aes(x=blob_end_D-end_se, y=y_fixed, xend=blob_end_D+end_se, yend=y_fixed),color="black")
)


# load latency data
dxy <- read.table("./data/dxy_all_rtbin.txt",sep=",",header=T)
d <- aggregate(rt~trial_id+bin+id,dxy,mean)
rm(dxy)

xdensity <-   ggplot(d, aes(x=-rt,group=id,fill=id))+ geom_vline(xintercept = 0, size = 0.3) + geom_density(position="stack",color="white",size=0,alpha=1)+ nice_theme+coord_cartesian(xlim=c(-950,-38),ylim=c(0,0.06))+ labs(x = " ",y="density")+scale_x_continuous(breaks=seq(-1000,0,200),labels=seq(-1000,0,200)/1000)+scale_y_continuous(breaks=c(0,0.1),labels=c("   0.0","0.100"))+facet_grid(.~bin)+scale_fill_grey(guide=F, start = 0, end=0.6)

plot_group_output <- grid.arrange(xdensity,  kernel_plot, ncol=1, nrow=2, heights=c(0.75,1.25))

ggsave("./figs/bayesian_rt_split_veraged_4bin.pdf", plot=plot_group_output, device="pdf", width=9,height=4)

