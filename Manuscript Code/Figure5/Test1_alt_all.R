#####################################
# Pooled radius
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


Calphas <- parsePDB("2i0q.pdb", subunit="A")
dssp <- read.dssp("2I0Q_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# Null scenario
#
# pooled 6.5, 9, 11.5 radius

set.seed(572019)
# Generate null

perm.reps <- 1000

reps<-1000 # ultimately want 1000
mean.p.65<-rep(NA,reps)
sd.p.65 <- rep(NA,reps)
ks.p.65 <- rep(NA,reps)
tail.p.65 <- rep(NA, reps)
quantile.p.65.all <- rep(NA, reps)
quantile.p.65.subs <- rep(NA, reps)

mean.p.9<-rep(NA,reps)
sd.p.9 <- rep(NA,reps)
ks.p.9 <- rep(NA,reps)
tail.p.9 <- rep(NA, reps)
quantile.p.9.all <- rep(NA, reps)
quantile.p.9.subs <- rep(NA, reps)

mean.p.115<-rep(NA,reps)
sd.p.115 <- rep(NA,reps)
ks.p.115 <- rep(NA,reps)
tail.p.115 <- rep(NA, reps)
quantile.p.115.all <- rep(NA, reps)
quantile.p.115.subs <- rep(NA, reps)

mean.p<-rep(NA,reps)
sd.p <- rep(NA,reps)
ks.p <- rep(NA,reps)
tail.p <- rep(NA, reps)
quantile.p.all <- rep(NA, reps)
quantile.p.subs <- rep(NA, reps)

yu.p <- rep(NA, reps)




for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)
  yu.p[l] <- yu.test(Calphas, subs, perm.reps)


  # do all versions of our test
  perm.subs.65 <- permute(Calphas, subs, reps=perm.reps, radius=6.5, contact=TRUE, similarity=1)
  perm.subs.9 <- permute(Calphas, subs, reps=perm.reps, radius=9, contact=TRUE, similarity=1)
  perm.subs.115 <- permute(Calphas, subs, reps=perm.reps, radius=11.5, contact=TRUE, similarity=1)
  
  perm.dist.65 <- apply(perm.subs.65, 2, fractions, Cdata = Calphas, radius = 6.5)
  perm.dist.9 <- apply(perm.subs.9, 2, fractions, Cdata = Calphas, radius = 9)
  perm.dist.115 <- apply(perm.subs.115, 2, fractions, Cdata = Calphas, radius = 11.5)

  perm.dist.65.subsonly <- apply(perm.subs.65, 2, fractions_subsonly, Cdata = Calphas, radius = 6.5)
  perm.dist.9.subsonly <- apply(perm.subs.9, 2, fractions_subsonly, Cdata = Calphas, radius = 9)
  perm.dist.115.subsonly <- apply(perm.subs.115, 2, fractions_subsonly, Cdata = Calphas, radius = 11.5)

  # radius 6.5
  perm.all <- as.vector(perm.dist.65)
  perm.mean <- apply(perm.dist.65, 2, mean)
  perm.sd <- apply(perm.dist.65, 2, sd)
  perm.95th.all <- apply(perm.dist.65, 2, quantile, probs = 0.95, type=2)
  perm.95th.subs <- apply(perm.dist.65.subsonly, 2, quantile, probs = 0.95, type=2)

  data.dist.65<-fractions(Calphas,subs,6.5)
  data.mean <- mean(data.dist.65)
  data.sd <- sd(data.dist.65)
  data.95th.all <- quantile(data.dist.65, 0.95, type=2)

  mean.p.65[l] <- sum(perm.mean > data.mean)/length(perm.mean)
  sd.p.65[l] <- sum(perm.sd > data.sd)/length(perm.sd)
  ks.p.65[l] <- ks.test(data.dist.65, perm.all, alternative="less")$p.value
  tail.p.65[l] <- tail.prop.test(data.dist.65, perm.all, 0.05)
  quantile.p.65.all[l] <- sum(perm.95th.all >= data.95th.all)/length(perm.95th.all)

  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.95th.subs <- quantile(data.dist.65.subsonly, 0.95)
  quantile.p.65.subs[l] <- sum(perm.95th.subs >= data.95th.subs)/length(perm.95th.subs)



  # radius 9
  perm.all <- as.vector(perm.dist.9)
  perm.mean <- apply(perm.dist.9, 2, mean)
  perm.sd <- apply(perm.dist.9, 2, sd)
  perm.95th.all <- apply(perm.dist.9, 2, quantile, probs = 0.95, type=2)
  perm.95th.subs <- apply(perm.dist.9.subsonly, 2, quantile, probs = 0.95, type=2)

  data.dist.9<-fractions(Calphas,subs,9)
  data.mean <- mean(data.dist.9)
  data.sd <- sd(data.dist.9)
  data.95th.all <- quantile(data.dist.9, 0.95, type=2)

  mean.p.9[l] <- sum(perm.mean > data.mean)/length(perm.mean)
  sd.p.9[l] <- sum(perm.sd > data.sd)/length(perm.sd)
  ks.p.9[l] <- ks.test(data.dist.9, perm.all, alternative="less")$p.value
  tail.p.9[l] <- tail.prop.test(data.dist.9, perm.all, 0.05)
  quantile.p.9.all[l] <- sum(perm.95th.all >= data.95th.all)/length(perm.95th.all)

  data.dist.9.subsonly <- fractions_subsonly(Calphas, subs, 9)
  data.95th.subs <- quantile(data.dist.9.subsonly, 0.95)
  quantile.p.9.subs[l] <- sum(perm.95th.subs >= data.95th.subs)/length(perm.95th.subs)


  # radius 11.5
  perm.all <- as.vector(perm.dist.115)
  perm.mean <- apply(perm.dist.115, 2, mean)
  perm.sd <- apply(perm.dist.115, 2, sd)
  perm.95th.all <- apply(perm.dist.115, 2, quantile, probs = 0.95, type=2)
  perm.95th.subs <- apply(perm.dist.115.subsonly, 2, quantile, probs = 0.95, type=2)

  data.dist.115<-fractions(Calphas,subs,6.5)
  data.mean <- mean(data.dist.115)
  data.sd <- sd(data.dist.115)
  data.95th.all <- quantile(data.dist.115, 0.95, type=2)

  mean.p.115[l] <- sum(perm.mean > data.mean)/length(perm.mean)
  sd.p.115[l] <- sum(perm.sd > data.sd)/length(perm.sd)
  ks.p.115[l] <- ks.test(data.dist.115, perm.all, alternative="less")$p.value
  tail.p.115[l] <- tail.prop.test(data.dist.115, perm.all, 0.05)
  quantile.p.115.all[l] <- sum(perm.95th.all >= data.95th.all)/length(perm.95th.all)

  data.dist.115.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.95th.subs <- quantile(data.dist.115.subsonly, 0.95)
  quantile.p.115.subs[l] <- sum(perm.95th.subs >= data.95th.subs)/length(perm.95th.subs)





  # pooled 
  perm.all <- c(as.vector(perm.dist.65),as.vector(perm.dist.9),as.vector(perm.dist.115))
  perm.all.stacked <- rbind(perm.dist.65, perm.dist.9, perm.dist.115)
  

  perm.mean <- apply(perm.all.stacked, 2, mean)
  perm.sd <- apply(perm.all.stacked, 2, sd)
  perm.95th.all <- apply(perm.all.stacked, 2, quantile, probs = 0.95, type=2)

  perm.subs.stacked <- rbind(perm.dist.65.subsonly, perm.dist.9.subsonly, perm.dist.115.subsonly)
  perm.95th.subs <- apply(perm.subs.stacked, 2, quantile, probs = 0.95, type=2)


  data.dist.65<-fractions(Calphas,subs,6.5)
  data.dist.9<-fractions(Calphas,subs,9)
  data.dist.115<-fractions(Calphas,subs,11.5)
  data.dist.pooled <- c(data.dist.65, data.dist.9, data.dist.115)

  data.mean <- mean(data.dist.pooled)
  data.sd <- sd(data.dist.pooled)
  data.95th.all <- quantile(data.dist.pooled, 0.95, type=2)

  mean.p[l] <- sum(perm.mean > data.mean)/length(perm.mean)
  sd.p[l] <- sum(perm.sd > data.sd)/length(perm.sd)
  ks.p[l] <- ks.test(data.dist.pooled, perm.all, alternative="less")$p.value
  tail.p[l] <- tail.prop.test(data.dist.pooled, perm.all, 0.05)
  quantile.p.all[l] <- sum(perm.95th.all >= data.95th.all)/length(perm.95th.all)

  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.dist.9.subsonly <- fractions_subsonly(Calphas, subs, 9)
  data.dist.115.subsonly <- fractions_subsonly(Calphas, subs, 11.5)
  data.subs.pooled <- c(data.dist.65.subsonly, data.dist.9.subsonly, data.dist.115.subsonly)

  data.95th.subs <- quantile(data.subs.pooled, 0.95)
  quantile.p.subs[l] <- sum(perm.95th.subs >= data.95th.subs)/length(perm.95th.subs)


}


save.image("test1_alt_all.RData")

