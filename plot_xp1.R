# --------------------------------------------------------------- #
# plot group level functions, exp. 1

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")

filelist <- dir("./results/xp1/")
ids <- {}; for(i in filelist) ids <- c(ids, strsplit(i,"_")[[1]][1])
ids <- unique(ids)

# plotting tools
library(ggplot2)

# this is only for loading a custom ggplot theme...
source("/home/matteo/sync/miscR/miscFunctions.R") 

# here is the theme I used:
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

library(mlisi)
library(rethinking)

# make group level dataset
group_dat <- {}
for(id_i in ids){
  file_saccade <- paste("./results/xp1/",id_i,"_saccade_xp1_fit.RDS",sep="")
  file_perception <- paste("./results/xp1/",id_i,"_perception_xp1_fit.RDS",sep="")
  m0 <- readRDS(file_saccade)
  m1 <- readRDS(file_perception)
  s_samples <- extract(m0, pars=c("beta"))$beta
  p_samples <- extract(m1, pars=c("beta"))$beta
  plo <- data.frame(w=c(apply(s_samples, 2, mean),apply(p_samples, 2, mean)),
                    w_se=c(apply(s_samples, 2, sd),apply(p_samples, 2, sd)),
                    lb = c(apply(s_samples, 2, function(x) HPDI(x, prob=0.95)[1]), apply(p_samples, 2, function(x) HPDI(x, prob=0.95)[1])),
                    ub = c(apply(s_samples, 2, function(x) HPDI(x, prob=0.95)[2]), apply(p_samples, 2, function(x) HPDI(x, prob=0.95)[2])),
                    D = rep(seq(-900,-9,length.out=100)+9/2,2),
                    type=c(rep("saccade",100),rep("perception",100)),
                    id = id_i)
  group_dat <- rbind(group_dat, plo)
}
str(group_dat)
N <- length(unique(group_dat$id))

# average and compute standard error across participants
dag <- aggregate(w~D+type,group_dat, mean)
dag$se <- aggregate(w~D+type,group_dat, se)$w
str(dag)

# set a probability threshold
alpha <- 0.05

# thresholding
group_dat$signif <- ifelse(group_dat$lb>0 | group_dat$ub<0,1,NA)

#cluster test - custom functions
source("cluster_test_functions.R")

## set parameters for cluster test
n_resels <- 900 / (1/15 * 1000)
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
    
    blob_p <- cluster_p(th, k=blob_length, S=100, R=n_resels)
    
    d_line <- data.frame(n=n_blob, s=blob_length, type=group_dat$type[i-1],p = blob_p, blob_start=start_blob, blob_end = end_blob, blob_start_D= start_blob_D, blob_end_D = end_blob_D, id=group_dat$id[i-1])
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


# change based on ID the 'height' of signi_plot horizontal lines
# then add "significance" column to the aggregate database
group_dat$id_n <- as.numeric(group_dat$id)
group_dat$signi_plot <- group_dat$signi* (-0.0135) - 0.004*group_dat$id_n
group_dat$se <- NA

# calculate the mean start and end points of the intergration windows
saccade_avg <- c(mean(d_blob_out$blob_start_D[d_blob_out$type=="saccade"]),mean(d_blob_out$blob_end_D[d_blob_out$type=="saccade"])) /1000
saccade_avg_se <- c(se(d_blob_out$blob_start_D[d_blob_out$type=="saccade"]),se(d_blob_out$blob_end_D[d_blob_out$type=="saccade"])) /1000

# for the perception I take and align only the bigger blob
blob_per <- d_blob_out[d_blob_out$type=="perception",]
perception_avg <- c(mean(blob_per$blob_start_D),mean(blob_per$blob_end_D)) /1000
perception_avg_se <- c(se(blob_per$blob_start_D),se(blob_per$blob_end_D)) /1000

# group y coordinate for segments
y_group <- (-0.02) - 0.004*(max(group_dat$id_n)+1)

group_plot <- ggplot(dag,aes(x=D/1000,y=w,ymin=w-se,ymax=w+se,color=type, fill=type))+geom_line(data=group_dat,aes(group=id:type),alpha=0.2)+geom_ribbon(alpha=0.25,color=NA)+geom_line(size=1.2)+geom_hline(yintercept=0,lty=2)+coord_cartesian(ylim=c(-0.045,0.11),xlim=c(-800,-38)/1000)+geom_vline(xintercept=0,lty=1)+scale_x_continuous(breaks=seq(-1000,0,100)/1000)+scale_y_continuous(breaks=seq(0,1,0.02))+nice_theme+scale_color_manual(values=c("black","blue"),name="")+scale_fill_manual(values=c("black","blue"),name="")+geom_line(data=group_dat, aes(y=signi_plot,group=id:type),size=1.2)+labs(x="time relative to saccade onset [sec]",y="weight")+annotate('segment',x=perception_avg[1],xend=perception_avg[2],y=y_group,yend=y_group,size=5, color="black",alpha=0.3)+annotate('segment',x=saccade_avg[1],xend=saccade_avg[2],y=y_group,yend=y_group,size=5, color="blue", alpha=0.3)+geom_segment(x=saccade_avg[1]-saccade_avg_se[1],xend=saccade_avg[1]+saccade_avg_se[1],y=y_group,yend=y_group,color="blue")+geom_segment(x=saccade_avg[2]-saccade_avg_se[2],xend=saccade_avg[2]+saccade_avg_se[2],y=y_group,yend=y_group, color="blue")+geom_segment(x=perception_avg[1]-perception_avg_se[1],xend=perception_avg[1]+perception_avg_se[1],y=y_group,yend=y_group, color="black")+geom_segment(x=perception_avg[2]-perception_avg_se[2],xend=perception_avg[2]+perception_avg_se[2],y=y_group,yend=y_group, color="black")

ggsave("./figs/group_xp1.pdf", plot=group_plot, device="pdf", width=6,height=2.5)




