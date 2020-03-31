# --------------------------------------------------------------- #
# "build" main dataset by applying preliminary operations such as cut-offs, etc.
# Matteo Lisi, 2019
#

rm(list=ls())
setwd("~/git_local/serial-integration-analysis/")

# --------------------------------------------------------------- #
# cut-offs
min_latency <- 100
min_amplitude <- 2.5

# --------------------------------------------------------------- #
# prepare data exp 1 (cut-offs, etc.)
d <- read.table("./data/raw_data_exp1.txt",header=T,sep="\t")
d <- d[d$sacType==1,] # keep only main response saccade
d$id <- factor(substr(d$vpcode,3,4))
d <- d[!is.na(d$id),] # some lines are unlabelled? (only 40 out of 140200); remove them

# get saccade components relative to starting location
d$sacXresp <- abs(d$sacxOffset - d$sacxOnset)
d$sacYresp <- abs(d$sacyOffset - d$sacyOnset)

dag <- aggregate(cbind(sacOnset, correct, errorChange, sacVPeak,sacXresp)~trial3+session+block+id, d, mean)
tapply(dag$id,dag$id,length) # observers run varying number of trials
# js   lo   ml   pc
# 988 1857 3219  944

mean(dag$sacOnset < min_latency) *100 # 0.4851598
mean(dag$sacXresp < min_amplitude) *100 # 5.065639
mean(dag$errorChange >0) *100 # 2.354452

# apply cut-offs
d <- d[d$sacOnset >= min_latency, ]
d <- d[d$errorChange<0,]
d <- d[d$sacXresp>=min_amplitude,]

# coordinates of choosen target
d$tarX <- ifelse(d$correct==1, d$meanXnear, d$meanXfar) # selected target
d$xs <- ifelse(d$correct==1, d$xNear, d$xFar) # xS indicate selected (saccaded to) target
d$xNs <- ifelse(d$correct==0, d$xNear, d$xFar)
d$ys <- ifelse(d$correct==1, d$yNear, d$yFar)

# compute saccadic error
d$errX <- d$sacXresp - abs(d$tarX)
d$errY <- d$sacyOffset

# response (R/L)
# side indicate the side of the near target (-1 = left; +1 = right)
d$resp <- ifelse((d$correct==1 & d$side==1)|(d$correct==0 & d$side==-1), 1, 0)

# transform coordinates for model fitting (disregard direction)
d$tarX <- abs(d$tarX)
d$xs <- abs(d$xs) - d$tarX # absolute values of xS (greater values indicate greater eccentricities)

# comput time of the stimulus samples relative to saccade onset
d$tSteps <- d$tSteps - d$sacOnset

write.table(d, file="./data/exp1_data",quote=F,row.names=F,sep="\t")

# --------------------------------------------------------------- #
# do the same for the data of Experiment 2
rm(d)
d <- read.table("./data/raw_data_exp2.txt",header=T,sep="\t")
d <- d[d$sacType==1,] # keep only main response saccade
d$id <- substr(d$vpcode,3,4)
d <- d[!is.na(d$id),] # some lines are unlabelled? (only 20 out of 80540); remove them
d$sacXresp <- abs(d$sacxOffset - d$sacxOnset)
d$sacYresp <- abs(d$sacyOffset - d$sacyOnset)
dag <- aggregate(cbind(sacOnset, correct, errorChange, sacVPeak,sacXresp)~trial3+session+block+id, d, mean)
tapply(dag$id,dag$id,length)
tapply(dag$correct,dag$id,mean)
tapply(dag$sacXresp,dag$id,mean)

mean(dag$sacOnset < min_latency) *100 # 0.3477397
mean(dag$sacXresp < min_amplitude) *100 # 20.59116
mean(dag$errorChange >0) *100 # 3.42772

d <- d[d$sacOnset >= min_latency, ]
d <- d[d$errorChange<0,]
d <- d[d$sacXresp>=min_amplitude,]

d$tarX <- ifelse(d$correct==1, d$meanXbright, d$meanXbase)
d$xs <- ifelse(d$correct==1, d$xBright, d$xBase)
d$ys <- ifelse(d$correct==1, d$yBright, d$yBase)
d$errX <- d$sacXresp - abs(d$tarX)
d$errY <- d$sacyOffset

# response (R/L)
d$resp <- ifelse((d$correct==1 & d$side==1)|(d$correct==0 & d$side==-1), 1, 0)

# transform for fitting
d$tarX <- abs(d$tarX)
d$xs <- abs(d$xs) - d$tarX # absolute values of xS !

# comput time of the samples relative to saccade onset
d$tSteps <- d$tSteps - d$sacOnset

# correct for side (such that differences are right minus left side)
d$LumDiff <- -d$side * d$LumDiff
d$meanLumDiff <- -d$side * d$meanLumDiff

# put on a common scale similar to the other experiment
d$LumDiff <- scale(d$LumDiff, center=F) * 2.5
d$meanLumDiff <- scale(d$meanLumDiff, center=F) * 2.5

# decorrelated (for fitting) istantaneous and stable luminance differences within a sequence
d$LumDiff <- d$LumDiff-d$meanLumDiff

write.table(d, file="./data/exp2_data",quote=F,row.names=F,sep="\t")

