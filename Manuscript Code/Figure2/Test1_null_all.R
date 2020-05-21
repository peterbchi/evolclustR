#####################################
# Radius 9
#
# 7/5/19 in progress

# setwd("C:/Users/pchi01/Dropbox/Temple/Clustering/analyses/Manuscript/April2019")
setwd("/home/peter/ClusterSignificance/Manuscript/Test1")

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
source("/home/peter/ClusterSignificance/evolClustR/R/permuteone.R")
source("/home/peter/ClusterSignificance/evolClustR/R/permute.R")
source("/home/peter/ClusterSignificance/evolClustR/R/permutetest.R")
source("/home/peter/ClusterSignificance/evolClustR/R/yutest.R")
source("/home/peter/ClusterSignificance/evolClustR/R/countsubs.R")


Calphas <- parsePDB("2i0q.pdb", subunit="A")
dssp <- read.dssp("2I0Q_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# Null scenario
#
# radius: 6.5, 9, 11
# and yu test

set.seed(572019)
# Generate null

# all this needs to be changed/fixed, still many remnants from parametric bootstrap
perm.reps <- 1000

reps<-1000 # ultimately want 1000
quantile.p.65.subs <- rep(NA, reps)
quantile.p.9.subs <- rep(NA, reps)
quantile.p.115.subs <- rep(NA, reps)
yu.p <- rep(NA, reps)

for(l in 1:reps){
  subs <- sample(Calphas$pos, 20)
  # perm.subs <- permute(Calphas, subs, reps=perm.reps, radius=9, contact=TRUE, similarity=1)
  # I don't know why this was here
  
  quantile.p.65.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=6.5, contact=FALSE, test=7)
  quantile.p.9.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=9, contact=FALSE, test=7)
  quantile.p.115.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=11.5, contact=FALSE, test=7)
  yu.p[l] <- yu.test(Calphas, subs, perm.reps)

  print(l)
save.image("/home/peter/ClusterSignificance/Manuscript/FinalCode/Figure2_mean/test1_null_all.RData")
}


save.image("/home/peter/ClusterSignificance/Manuscript/FinalCode/Figure2_mean/test1_null_all.RData")



