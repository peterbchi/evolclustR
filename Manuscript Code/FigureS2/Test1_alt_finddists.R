#####################################
# 
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
# 2i0q
#

set.seed(572019)
# Generate null
reps <- 1000

dists_2i0q <- matrix(NA,nrow=reps,ncol=190)
for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)
  Cdata_sub <- Calphas[is.element(Calphas$pos, subs),]
  tmp <- sapply(1:20, find.dist, Cdata=Cdata_sub)
  dists_2i0q[l,] <- tmp[upper.tri(tmp)]
  print(l)
}
dists_2i0q <- as.vector(dists_2i0q)


####################
# 1D4t

setwd("...")

Calphas <- parsePDB("1d4t.pdb", subunit="A")
dssp <- read.dssp("1D4T_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

dists_1d4t <- matrix(NA,nrow=reps,ncol=10)
for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 5, ratio = 1, exponent = 2)
  Cdata_sub <- Calphas[is.element(Calphas$pos, subs),]
  tmp <- sapply(1:5, find.dist, Cdata=Cdata_sub)
  dists_1d4t[l,] <- tmp[upper.tri(tmp)]
  print(l)
}
dists_1d4t <- as.vector(dists_1d4t)


################
# 1ax8

setwd("...")
Calphas <- parsePDB("1ax8.pdb", subunit="A")
dssp <- read.dssp("1ax8_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)



dists_1ax8 <- matrix(NA,nrow=reps,ncol=190)
for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 20, ratio = 1, exponent = 2)
  Cdata_sub <- Calphas[is.element(Calphas$pos, subs),]
  tmp <- sapply(1:20, find.dist, Cdata=Cdata_sub)
  dists_1ax8[l,] <- tmp[upper.tri(tmp)]
  print(l)
}
dists_1ax8 <- as.vector(dists_1ax8)



setwd("...")
save.image("test1_alt_all_dists.RData")

