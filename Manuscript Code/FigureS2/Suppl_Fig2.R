setwd("...")


load("test1_alt_all_dists.RData")

pdf(".../Suppl_Fig2.pdf")
par(mfrow=c(3,1))
hist(dists_2i0q,main="2I0Q",xlim=c(0,100),nclass=100,xlab="")
hist(dists_1d4t,main="1D4T",xlim=c(0,100),nclass=100,xlab="")
hist(dists_1ax8,main="1AX8",xlim=c(0,100),nclass=100,xlab="Angstroms")
dev.off()



