#####################################
# 3 on surface
# Following G118 from Gaucher2003
#

# setwd("C:/Users/pchi01/Dropbox/Temple/Clustering/analyses/Manuscript/April2019")
setwd("/home/peter/ClusterSignificance/Manuscript/Leptin")

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
source("/home/peter/ClusterSignificance/evolClustR/R/permutetest.R")
source("/home/peter/ClusterSignificance/evolClustR/R/readsubrecon.R")


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

boot.95th.subs.65 <- rep(NA, boot.reps)


for(i in 1:boot.reps){
  boot.subs <- generate.sasa.subs(Calphas, b.length = 6, ratio = 1, exponent = 2)
  boot.dist.65.subsonly <- fractions_subsonly(Calphas, boot.subs, 6.5)
  boot.95th.subs.65[i] <- quantile(boot.dist.65.subsonly, 0.95)  

  print(i)
}


reps<-1000 # ultimately want 1000

quantile.p.65.subs <- rep(NA, reps)
yu.p <- rep(NA, reps)

dists_1ax8 <- matrix(NA,nrow=reps,ncol=15)


for(l in 1:reps){
  # initialize three deterministic substitutions
  # first one: G118
  subs <- rep(0, dim(Calphas)[1])
  most.acc <- which(Calphas$pos == 118)
  subs[most.acc] <- 1

  # second one: closest to first one
  dists <- find.dist(Calphas, most.acc)
  dists[which(dists==0)] <- max(dists)
  closest <- which(dists == min(dists))
  subs[closest] <- 1
  
  # third one: next closest
  dists[closest] <- max(dists)
  next.closest <- which(dists == min(dists))
  subs[next.closest] <- 1
  

  subs <- generate.sasa.subs(Calphas, b.length = 6, ratio = 1, exponent = 2, subs = subs)

  # do Yu/Thorne test
  yu.p[l] <- yu.test(Calphas, subs, boot.reps)

  # store pairwise distances between subs (still fix)
  Cdata_sub <- Calphas[is.element(Calphas$pos, subs),]
  tmp <- sapply(1:6, find.dist, Cdata=Cdata_sub)
  dists_1ax8[l,] <- tmp[upper.tri(tmp)]
  
  
  data.dist.65.subsonly <- fractions_subsonly(Calphas, subs, 6.5)
  data.95th.subs.65 <- quantile(data.dist.65.subsonly, 0.95)  

  quantile.p.65.subs[l] <- sum(boot.95th.subs.65 >= data.95th.subs.65)/length(boot.95th.subs.65)

  print(l)
  save.image("test2_surface3_macaca_1ax8.RData")
}

macaca_subs <- c(8, 48, 53, 66, 67, 112)
Cdata_sub <- Calphas[is.element(Calphas$pos, macaca_subs),]
tmp <- sapply(1:6, find.dist, Cdata=Cdata_sub)
dists_hominoid_1ax8 <- tmp[upper.tri(tmp)]



save.image("test2_surface3_macaca_1ax8.RData")

