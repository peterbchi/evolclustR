#####################################
# 2 on surface
#
#

# setwd("C:/Users/pchi01/Dropbox/Temple/Clustering/analyses/Manuscript/April2019")
setwd("/home/peter/ClusterSignificance/Manuscript/Test2")

source("/home/peter/ClusterSignificance/evolClustR/R/countneighbors.R")
source("/home/peter/ClusterSignificance/evolClustR/R/fractions.R")
source("/home/peter/ClusterSignificance/evolClustR/R/findneighbors.R")
source("/home/peter/ClusterSignificance/evolClustR/R/finddist.R")
source("/home/peter/ClusterSignificance/evolClustR/R/readdssp.R")
source("/home/peter/ClusterSignificance/evolClustR/R/generatesasasubs.R")
source("/home/peter/ClusterSignificance/evolClustR/R/parsePDB.R")
source("/home/peter/ClusterSignificance/evolClustR/R/calcsasa.R")
source("/home/peter/ClusterSignificance/evolClustR/R/tailproptest.R")
source("/home/peter/ClusterSignificance/evolClustR/R/fractions_subsonly.R")
source("/home/peter/ClusterSignificance/evolClustR/R/countsubs.R")
source("/home/peter/ClusterSignificance/evolClustR/R/yutest.R")
source("/home/peter/ClusterSignificance/evolClustR/R/permuteone.R")
source("/home/peter/ClusterSignificance/evolClustR/R/permute.R")


Calphas <- parsePDB("2i0q.pdb", subunit="A")
dssp <- read.dssp("2I0Q_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# 2 surface
# 6.5 to 9 radius
# 

set.seed(1032019)
# Generate null
boot.reps <- 1000

boot.95th.subs.65 <- rep(NA, boot.reps)
boot.95th.subs.7 <- rep(NA, boot.reps)
boot.95th.subs.75 <- rep(NA, boot.reps)
boot.95th.subs.8 <- rep(NA, boot.reps)
boot.95th.subs.85 <- rep(NA, boot.reps)
boot.95th.subs.9 <- rep(NA, boot.reps)


for(i in 1:boot.reps){
  boot.subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)

  boot.dist.65.subsonly <- fractions_subsonly(Calphas, boot.subs, 6.5)
  boot.dist.7.subsonly <- fractions_subsonly(Calphas, boot.subs, 7)
  boot.dist.75.subsonly <- fractions_subsonly(Calphas, boot.subs, 7.5)
  boot.dist.8.subsonly <- fractions_subsonly(Calphas, boot.subs, 8)
  boot.dist.85.subsonly <- fractions_subsonly(Calphas, boot.subs, 8.5)
  boot.dist.9.subsonly <- fractions_subsonly(Calphas, boot.subs, 9)

  boot.95th.subs.65[i] <- quantile(boot.dist.65.subsonly, 0.95)  
  boot.95th.subs.7[i] <- quantile(boot.dist.7.subsonly, 0.95)  
  boot.95th.subs.75[i] <- quantile(boot.dist.75.subsonly, 0.95)  
  boot.95th.subs.8[i] <- quantile(boot.dist.8.subsonly, 0.95)  
  boot.95th.subs.85[i] <- quantile(boot.dist.85.subsonly, 0.95)  
  boot.95th.subs.9[i] <- quantile(boot.dist.9.subsonly, 0.95)  



  print(i)
  save.image("test2_surface2_rads.RData")
}

load("test2_surface2_rads.RData")

reps<-1000 # ultimately want 1000

quantile.p.65.subs <- rep(NA, reps)
quantile.p.7.subs <- rep(NA, reps)
quantile.p.75.subs <- rep(NA, reps)
quantile.p.8.subs <- rep(NA, reps)
quantile.p.85.subs <- rep(NA, reps)
quantile.p.9.subs <- rep(NA, reps)

yu.p <- rep(NA, reps)



for(l in 1:reps){
  # initialize two deterministic substitutions
  # first one: highest acc.ratio as before
  subs <- rep(0, dim(Calphas)[1])
  most.acc <- which(Calphas$acc.ratio == max(Calphas$acc.ratio))
  subs[most.acc] <- 1

  # second one: closest to first one
  dists <- find.dist(Calphas, most.acc)
  dists[which(dists==0)] <- max(dists)
  closest <- which(dists == min(dists))
  subs[closest] <- 1

  subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2, subs = subs)

  # do Yu/Thorne test
  yu.p[l] <- yu.test(Calphas, subs, boot.reps)


  
  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.dist.7.subsonly <- fractions_subsonly(Calphas, subs, 7)
  data.dist.75.subsonly <- fractions_subsonly(Calphas, subs, 7.5)
  data.dist.8.subsonly <- fractions_subsonly(Calphas, subs, 8)
  data.dist.85.subsonly <- fractions_subsonly(Calphas, subs, 8.5)
  data.dist.9.subsonly <- fractions_subsonly(Calphas, subs, 9)

  data.95th.subs.65 <- quantile(data.dist.65.subsonly, 0.95)  
  data.95th.subs.7 <- quantile(data.dist.7.subsonly, 0.95)  
  data.95th.subs.75 <- quantile(data.dist.75.subsonly, 0.95)  
  data.95th.subs.8 <- quantile(data.dist.8.subsonly, 0.95)  
  data.95th.subs.85 <- quantile(data.dist.85.subsonly, 0.95)  
  data.95th.subs.9 <- quantile(data.dist.9.subsonly, 0.95)  


  quantile.p.65.subs[l] <- sum(boot.95th.subs.65 >= data.95th.subs.65)/length(boot.95th.subs.65)
  quantile.p.7.subs[l] <- sum(boot.95th.subs.7 >= data.95th.subs.7)/length(boot.95th.subs.7)
  quantile.p.75.subs[l] <- sum(boot.95th.subs.75 >= data.95th.subs.75)/length(boot.95th.subs.75)
  quantile.p.8.subs[l] <- sum(boot.95th.subs.8 >= data.95th.subs.8)/length(boot.95th.subs.8)
  quantile.p.85.subs[l] <- sum(boot.95th.subs.85 >= data.95th.subs.85)/length(boot.95th.subs.85)
  quantile.p.9.subs[l] <- sum(boot.95th.subs.9 >= data.95th.subs.9)/length(boot.95th.subs.9)


  print(l)
  save.image("test2_surface2_rads.RData")
}


save.image("test2_surface2_rads.RData")


