# --------------------------------------------------------------- #
# plot individual functions, - both exp

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")

filelist <- dir("./results/xp1/")
ids <- {}; for(i in filelist) ids <- c(ids, strsplit(i,"_")[[1]][1])
ids <- unique(ids)

# plotting tools
library(ggplot2)
library(see)
library(ggpubr)

# this is only for loading a custom ggplot theme...
source("/home/matteo/sync/miscR/miscFunctions.R") 

# here is the theme I used:
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))


# --------------------------------------------------------------- #
# first plot individual functions for EXP1
plot_individual <- function(id, store=1){
  
  file_saccade <- paste("./results/xp1/",id,"_saccade_xp1_fit.RDS",sep="")
  file_perception <- paste("./results/xp1/",id,"_perception_xp1_fit.RDS",sep="")
  
  plot_name <- paste("./figs/xp1/",id,"_xp1.pdf",sep="")
  plot_title <- paste("Exp. 1,",id)
  
  m0 <- readRDS(file_saccade)
  m1 <- readRDS(file_perception)
  
  require(rethinking)
  
  beta_samples <- extract(m0, pars=c("beta"))$beta
  plo <- data.frame(beta=apply(beta_samples, 2, mean))
  plo$lb <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[1])
  plo$ub <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[2])
  plo$lb_se <- plo$beta - apply(beta_samples, 2, sd)
  plo$ub_se <- plo$beta + apply(beta_samples, 2, sd)
  plo$D <- seq(-1000,-10,10)+5
  plo$signif <- ifelse(plo$lb>0 | plo$ub<0,1,0)
  plo$signif <- ifelse(plo$signif, plo$signif * (-0.04), NA)
  plo$type <- "saccade"
  
  beta2_samples <- extract(m1, pars=c("beta"))$beta
  plo2 <- data.frame(beta=apply(beta2_samples, 2, mean))
  plo2$lb <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[1])
  plo2$ub <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[2])
  plo2$lb_se <- plo2$beta - apply(beta2_samples, 2, sd)
  plo2$ub_se <- plo2$beta + apply(beta2_samples, 2, sd)
  plo2$D <- seq(-1000,-10,10)+5
  plo2$signif <- ifelse(plo2$lb>0 | plo2$ub<0,1,0)
  plo2$signif <- ifelse(plo2$signif, plo2$signif * (-0.04), NA)
  plo2$type <- "perception"
  
  plo <- rbind(plo,plo2)
  
  outplot <- ggplot(plo,aes(x=D/1000,y=beta,ymin=lb,ymax=ub,color=type, fill=type))+geom_ribbon(alpha=0.25,color=NA)+geom_line(size=1.2)+geom_hline(yintercept=0,lty=2)+coord_cartesian(ylim=c(-0.05,0.20),xlim=c(-800,-38)/1000)+geom_vline(xintercept=0,lty=1)+scale_x_continuous(breaks=seq(-1000,0,100)/1000)+nice_theme+scale_color_manual(values=c("black","blue"),name="",guide=F)+scale_fill_manual(values=c("black","blue"),name="",guide=F)+geom_line(aes(y=signif),size=2)+labs(x="time relative to saccade onset [sec]",y="weight")+ggtitle(plot_title)
  
  if(store==1){
    ggsave(plot_name, plot=outplot, device="pdf", width=4,height=2)
  }
  print(outplot)
}

a <- plot_individual(ids[1], store=F)
b <- plot_individual(ids[2], store=F)
c <- plot_individual(ids[3], store=F)
d <- plot_individual(ids[4], store=F)

# --------------------------------------------------------------- #
# then add those for exp 2

filelist <- dir("./results/xp2/")
ids <- {}; for(i in filelist) ids <- c(ids, strsplit(i,"_")[[1]][1])
ids <- unique(ids)

plot_individual <- function(id, store=1,display=0){
  
  file_saccade <- paste("./results/xp2/",id,"_saccade_xp2_fit.RDS",sep="")
  file_perception <- paste("./results/xp2/",id,"_perception_xp2_fit.RDS",sep="")
  
  plot_name <- paste("./figs/xp2/",id,"_xp2.pdf",sep="")
  plot_title <- paste("Exp. 2,",id)
  
  
  m0 <- readRDS(file_saccade)
  m1 <- readRDS(file_perception)
  
  require(rethinking)
  
  beta_samples <- extract(m0, pars=c("beta"))$beta
  plo <- data.frame(beta=apply(beta_samples, 2, mean))
  plo$lb <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[1])
  plo$ub <- apply(beta_samples, 2, function(x) HPDI(x, prob=0.95)[2])
  plo$lb_se <- plo$beta - apply(beta_samples, 2, sd)
  plo$ub_se <- plo$beta + apply(beta_samples, 2, sd)
  plo$D <- seq(-1000,-10,10)+5
  plo$signif <- ifelse(plo$lb>0 | plo$ub<0,1,0)
  plo$signif <- ifelse(plo$signif, plo$signif * (-0.04), NA)
  plo$type <- "saccade"
  
  beta2_samples <- extract(m1, pars=c("beta"))$beta
  plo2 <- data.frame(beta=apply(beta2_samples, 2, mean))
  plo2$lb <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[1])
  plo2$ub <- apply(beta2_samples, 2, function(x) HPDI(x, prob=0.95)[2])
  plo2$lb_se <- plo2$beta - apply(beta2_samples, 2, sd)
  plo2$ub_se <- plo2$beta + apply(beta2_samples, 2, sd)
  plo2$D <- seq(-1000,-10,10)+5
  plo2$signif <- ifelse(plo2$lb>0 | plo2$ub<0,1,0)
  plo2$signif <- ifelse(plo2$signif, plo2$signif * (-0.04), NA)
  plo2$type <- "perception"
  
  plo <- rbind(plo,plo2)
  
  outplot <- ggplot(plo,aes(x=D/1000,y=beta,ymin=lb,ymax=ub,color=type, fill=type))+geom_ribbon(alpha=0.25,color=NA)+geom_line(size=1.2)+geom_hline(yintercept=0,lty=2)+coord_cartesian(ylim=c(-0.05,0.20),xlim=c(-800,-38)/1000)+geom_vline(xintercept=0,lty=1)+scale_x_continuous(breaks=seq(-1000,0,100)/1000)+nice_theme+scale_color_manual(values=c("black","blue"),name="",guide=F)+scale_fill_manual(values=c("black","blue"),name="",guide=F)+geom_line(aes(y=signif),size=2)+labs(x="time relative to saccade onset [sec]",y="weight")+ggtitle(plot_title)
  
  if(store==1){
    ggsave(plot_name, plot=outplot, device="pdf", width=4,height=2)
  }
  if(display==1){
    print(outplot) 
  }
}

e <- plot_individual(ids[1], store=F, display=1)
f <- plot_individual(ids[2], store=F, display=1)
g <- plot_individual(ids[3], store=F, display=1)
h <- plot_individual(ids[4], store=F, display=1)
i <- plot_individual(ids[5], store=F, display=1)

# --------------------------------------------------------------- #
# make big plot
pdf("./figs/individual_all.pdf",width=12,height=6)
ggarrange(a,b,c,d, e, f, g, h, i,
          labels = c("A", "B","C","D","E","F","G","H","I"),
          #labels = ids,
          ncol = 3, nrow = 3)
dev.off()

