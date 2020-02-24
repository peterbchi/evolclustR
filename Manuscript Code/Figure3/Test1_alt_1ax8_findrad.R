#####################################
# Finding max power
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
source(".../evolClustR/R/permutetest.R")


Calphas <- parsePDB("1ax8.pdb", subunit="A")
dssp <- read.dssp("1ax8_dssp.txt")
Calphas <- calc.sasa(Calphas, dssp)

#############################
# 
set.seed(572019)
# Generate null

perm.reps <- 1000
reps <- 1000 # ultimately want 1000

yu.p <- rep(NA, reps)

radii <- seq(5, 20, by=0.5)
pvals <- matrix(NA, nrow=reps, ncol=length(radii))

for(l in 1:reps){
  subs <- generate.sasa.subs(Calphas, b.length = 6, ratio = 1, exponent = 2) # 130 total sites; 6/130=0.046
  yu.p[l] <- yu.test(Calphas, subs, perm.reps)

  for(r in 1:length(radii)){
    pvals[l,r] <- permute.test(Calphas, subs, reps=perm.reps, radius = radii[r], test=7)
  }

  print(l)
}


save.image(".../test1_alt_1ax8_findrad.RData")

