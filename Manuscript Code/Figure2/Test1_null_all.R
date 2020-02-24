#####################################
# Radius 9
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
source(".../evolClustR/R/permuteone.R")
source(".../evolClustR/R/permute.R")
source(".../evolClustR/R/permutetest.R")
source(".../evolClustR/R/yutest.R")
source(".../evolClustR/R/countsubs.R")


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

perm.reps <- 1000

reps<-1000 
quantile.p.65.subs <- rep(NA, reps)
quantile.p.9.subs <- rep(NA, reps)
quantile.p.115.subs <- rep(NA, reps)
yu.p <- rep(NA, reps)

for(l in 1:reps){
  subs <- sample(Calphas$pos, 20)

  quantile.p.65.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=6.5, contact=FALSE, test=6)
  quantile.p.9.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=9, contact=FALSE, test=6)
  quantile.p.115.subs[l] <- permute.test(Calphas, subs, reps=perm.reps, radius=11.5, contact=FALSE, test=6)
  yu.p[l] <- yu.test(Calphas, subs, perm.reps)

  print(l)
}


save.image(".../test1_null_all.RData")



