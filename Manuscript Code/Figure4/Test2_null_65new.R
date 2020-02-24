#####################################
# Radius 6.5
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


Calphas <- parsePDB("2i0q.pdb", subunit="A")
dssp <- read.dssp("2I0Q_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# Null scenario
#
# radius: 6.5

set.seed(572019)
# Generate null
boot.reps <- 1000
boot.mean <- rep(NA,boot.reps)
boot.sd <- rep(NA, boot.reps)
boot.95th.all <- rep(NA, boot.reps)
boot.95th.subs <- rep(NA, boot.reps)
boot.all <- NULL

for(i in 1:boot.reps){
  boot.subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)
  boot.dist.65 <- fractions(Calphas, boot.subs, 6.5)
  boot.all <- c(boot.all, boot.dist.65)
  boot.mean[i] <- mean(boot.dist.65)
  boot.sd[i] <- sd(boot.dist.65)
  boot.95th.all[i] <- quantile(boot.dist.65, 0.95)

  boot.dist.65.subsonly <- fractions_subsonly(Calphas, boot.subs, 6.5)
  boot.95th.subs[i] <- quantile(boot.dist.65.subsonly, 0.95)  

}



reps<-1000 # ultimately want 1000
mean.p.65<-rep(NA,reps)
sd.p.65 <- rep(NA,reps)
ks.p.65 <- rep(NA,reps)
tail.p.65 <- rep(NA, reps)
quantile.p.65.all <- rep(NA, reps)
quantile.p.65.subs <- rep(NA, reps)


for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)
  data.dist.65<-fractions(Calphas,subs,6.5)
  data.mean <- mean(data.dist.65)
  data.sd <- sd(data.dist.65)
  data.95th.all <- quantile(data.dist.65, 0.95)

  mean.p.65[l] <- sum(boot.mean > data.mean)/length(boot.mean)
  sd.p.65[l] <- sum(boot.sd > data.sd)/length(boot.sd)
  ks.p.65[l] <- ks.test(data.dist.65, boot.all, alternative="less")$p.value
  tail.p.65[l] <- tail.prop.test(data.dist.65, boot.all, 0.05)
  quantile.p.65.all[l] <- sum(boot.95th.all > data.95th.all)/length(boot.95th.all)

  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.95th.subs <- quantile(data.dist.65.subsonly, 0.95)
  quantile.p.65.subs[l] <- sum(boot.95th.subs > data.95th.subs)/length(boot.95th.subs)

}


save.image("test2_null_65new.RData")

