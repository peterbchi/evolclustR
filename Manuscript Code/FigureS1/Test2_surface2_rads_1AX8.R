#####################################
# 2 on surface
#
#

setwd("...")

source(".../evolClustR/R/countneighbors.R")
source(".../evolClustR/R/fractions.R")
source(".../evolClustR/R/findneighbors.R")
source(".../evolClustR/R/finddist.R")
source(".../evolClustR/R/readdssp.R")
source(".../evolClustR/R/generatesasasubs.R")
source(".../evolClustR/R/parsePDB.R")
source(".../evolClustR/R/calcsasa.R")
source(".../evolClustR/R/tailproptest.R")
source(".../evolClustR/R/fractions_subsonly.R")
source(".../evolClustR/R/countsubs.R")
source(".../evolClustR/R/yutest.R")
source(".../evolClustR/R/permuteone.R")
source(".../evolClustR/R/permute.R")


Calphas <- parsePDB("1ax8.pdb", subunit="A")
dssp <- read.dssp("1ax8_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# 4 surface
# 4.5 to 10 radius
# 

set.seed(1032019)
# Generate null
boot.reps <- 1000

boot.95th.subs.45 <- rep(NA, boot.reps)
boot.95th.subs.5 <- rep(NA, boot.reps)
boot.95th.subs.55 <- rep(NA, boot.reps)
boot.95th.subs.6 <- rep(NA, boot.reps)
boot.95th.subs.65 <- rep(NA, boot.reps)
boot.95th.subs.7 <- rep(NA, boot.reps)
boot.95th.subs.75 <- rep(NA, boot.reps)
boot.95th.subs.8 <- rep(NA, boot.reps)
boot.95th.subs.85 <- rep(NA, boot.reps)
boot.95th.subs.9 <- rep(NA, boot.reps)
boot.95th.subs.95 <- rep(NA, boot.reps)
boot.95th.subs.10 <- rep(NA, boot.reps)


for(i in 1:boot.reps){
  boot.subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)

  boot.dist.45.subsonly <- fractions_subsonly(Calphas, boot.subs, 4.5)
  boot.dist.5.subsonly <- fractions_subsonly(Calphas, boot.subs, 5)
  boot.dist.55.subsonly <- fractions_subsonly(Calphas, boot.subs, 5.5)
  boot.dist.6.subsonly <- fractions_subsonly(Calphas, boot.subs, 6)
  boot.dist.65.subsonly <- fractions_subsonly(Calphas, boot.subs, 6.5)
  boot.dist.7.subsonly <- fractions_subsonly(Calphas, boot.subs, 7)
  boot.dist.75.subsonly <- fractions_subsonly(Calphas, boot.subs, 7.5)
  boot.dist.8.subsonly <- fractions_subsonly(Calphas, boot.subs, 8)
  boot.dist.85.subsonly <- fractions_subsonly(Calphas, boot.subs, 8.5)
  boot.dist.9.subsonly <- fractions_subsonly(Calphas, boot.subs, 9)
  boot.dist.95.subsonly <- fractions_subsonly(Calphas, boot.subs, 9.5)
  boot.dist.10.subsonly <- fractions_subsonly(Calphas, boot.subs, 10)

  boot.95th.subs.45[i] <- quantile(boot.dist.45.subsonly, 0.95)  
  boot.95th.subs.5[i] <- quantile(boot.dist.5.subsonly, 0.95)  
  boot.95th.subs.55[i] <- quantile(boot.dist.55.subsonly, 0.95)  
  boot.95th.subs.6[i] <- quantile(boot.dist.6.subsonly, 0.95)  
  boot.95th.subs.65[i] <- quantile(boot.dist.65.subsonly, 0.95)  
  boot.95th.subs.7[i] <- quantile(boot.dist.7.subsonly, 0.95)  
  boot.95th.subs.75[i] <- quantile(boot.dist.75.subsonly, 0.95)  
  boot.95th.subs.8[i] <- quantile(boot.dist.8.subsonly, 0.95)  
  boot.95th.subs.85[i] <- quantile(boot.dist.85.subsonly, 0.95)  
  boot.95th.subs.9[i] <- quantile(boot.dist.9.subsonly, 0.95)  
  boot.95th.subs.95[i] <- quantile(boot.dist.95.subsonly, 0.95)  
  boot.95th.subs.10[i] <- quantile(boot.dist.10.subsonly, 0.95)  

}


reps<-1000 # ultimately want 1000

quantile.p.45.subs <- rep(NA, reps)
quantile.p.5.subs <- rep(NA, reps)
quantile.p.55.subs <- rep(NA, reps)
quantile.p.6.subs <- rep(NA, reps)
quantile.p.65.subs <- rep(NA, reps)
quantile.p.7.subs <- rep(NA, reps)
quantile.p.75.subs <- rep(NA, reps)
quantile.p.8.subs <- rep(NA, reps)
quantile.p.85.subs <- rep(NA, reps)
quantile.p.9.subs <- rep(NA, reps)
quantile.p.95.subs <- rep(NA, reps)
quantile.p.10.subs <- rep(NA, reps)

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


  
  data.dist.45.subsonly <- fractions_subsonly(Calphas, subs, 4.5)
  data.dist.5.subsonly <- fractions_subsonly(Calphas, subs, 5)
  data.dist.55.subsonly <- fractions_subsonly(Calphas, subs, 5.5)
  data.dist.6.subsonly <- fractions_subsonly(Calphas, subs, 6)
  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.dist.7.subsonly <- fractions_subsonly(Calphas, subs, 7)
  data.dist.75.subsonly <- fractions_subsonly(Calphas, subs, 7.5)
  data.dist.8.subsonly <- fractions_subsonly(Calphas, subs, 8)
  data.dist.85.subsonly <- fractions_subsonly(Calphas, subs, 8.5)
  data.dist.9.subsonly <- fractions_subsonly(Calphas, subs, 9)
  data.dist.95.subsonly <- fractions_subsonly(Calphas, subs, 9.5)
  data.dist.10.subsonly <- fractions_subsonly(Calphas, subs, 10)

  data.95th.subs.45 <- quantile(data.dist.45.subsonly, 0.95)  
  data.95th.subs.5 <- quantile(data.dist.5.subsonly, 0.95)  
  data.95th.subs.55 <- quantile(data.dist.55.subsonly, 0.95)  
  data.95th.subs.6 <- quantile(data.dist.6.subsonly, 0.95)  
  data.95th.subs.65 <- quantile(data.dist.65.subsonly, 0.95)  
  data.95th.subs.7 <- quantile(data.dist.7.subsonly, 0.95)  
  data.95th.subs.75 <- quantile(data.dist.75.subsonly, 0.95)  
  data.95th.subs.8 <- quantile(data.dist.8.subsonly, 0.95)  
  data.95th.subs.85 <- quantile(data.dist.85.subsonly, 0.95)  
  data.95th.subs.9 <- quantile(data.dist.9.subsonly, 0.95)  
  data.95th.subs.95 <- quantile(data.dist.95.subsonly, 0.95)  
  data.95th.subs.10 <- quantile(data.dist.10.subsonly, 0.95)  


  quantile.p.45.subs[l] <- sum(boot.95th.subs.45 >= data.95th.subs.45)/length(boot.95th.subs.45)
  quantile.p.5.subs[l] <- sum(boot.95th.subs.5 >= data.95th.subs.5)/length(boot.95th.subs.5)
  quantile.p.55.subs[l] <- sum(boot.95th.subs.55 >= data.95th.subs.55)/length(boot.95th.subs.55)
  quantile.p.6.subs[l] <- sum(boot.95th.subs.6 >= data.95th.subs.6)/length(boot.95th.subs.6)
  quantile.p.65.subs[l] <- sum(boot.95th.subs.65 >= data.95th.subs.65)/length(boot.95th.subs.65)
  quantile.p.7.subs[l] <- sum(boot.95th.subs.7 >= data.95th.subs.7)/length(boot.95th.subs.7)
  quantile.p.75.subs[l] <- sum(boot.95th.subs.75 >= data.95th.subs.75)/length(boot.95th.subs.75)
  quantile.p.8.subs[l] <- sum(boot.95th.subs.8 >= data.95th.subs.8)/length(boot.95th.subs.8)
  quantile.p.85.subs[l] <- sum(boot.95th.subs.85 >= data.95th.subs.85)/length(boot.95th.subs.85)
  quantile.p.9.subs[l] <- sum(boot.95th.subs.9 >= data.95th.subs.9)/length(boot.95th.subs.9)
  quantile.p.95.subs[l] <- sum(boot.95th.subs.95 >= data.95th.subs.95)/length(boot.95th.subs.95)
  quantile.p.10.subs[l] <- sum(boot.95th.subs.10 >= data.95th.subs.10)/length(boot.95th.subs.10)

}


save.image("test2_surface2_rads_1ax8.RData")

