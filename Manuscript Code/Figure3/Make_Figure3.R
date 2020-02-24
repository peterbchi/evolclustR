setwd("...")

pdf("Fig3_radius.pdf")
par(mfrow=c(3,1),oma=c(3,2.5,1.5,1),mar=c(3,2,1.5,1))

# 2I0Q
load("test1_alt_2i0q_findrad.RData")

power.temp <- pvals<0.05
power.2i0q <- apply(power.temp, 2, sum, na.rm=T)/sum(!is.na(power.temp[,1]))
plot(power.2i0q~radii, type='l',ylim=c(0,1), ylab="power", xlab="", main="2I0Q", cex.main=2, xaxt="n")
legend("topright", legend=c("CKL20","YT06"),col=c("black","gray"), lty=1, cex=1.5)
axis(1, at=5:20)

mtext("power",side=2,line=2.5,cex=1.2)

abline(h=sum(yu.p<0.05, na.rm=T)/sum(!is.na(yu.p)), col="gray")


# 1D4T
load("test1_alt_1d4t_findrad.RData")

power.temp <- pvals<0.05
power.1d4t <- apply(power.temp, 2, sum, na.rm=T)/sum(!is.na(power.temp[,1]))
plot(power.1d4t~radii, type='l',ylim=c(0,1), ylab="power", xlab="", main="1D4T", cex.main=2, xaxt="n")
axis(1, at=5:20)

mtext("power",side=2,line=2.5,cex=1.2)

abline(h=sum(yu.p<0.05)/length(yu.p), col="gray")


# 1AX8
load("test1_alt_1ax8_findrad.RData")

power.temp <- pvals<0.05
power.1ax8 <- apply(power.temp, 2, sum, na.rm=T)/sum(!is.na(power.temp[,1]))
plot(power.1ax8~radii, type='l',ylim=c(0,1), ylab="power", xlab="", main="1AX8", cex.main=2, xaxt="n")
axis(1, at=5:20)
mtext("radius",side=1,line=3,cex=1.2)
mtext("power",side=2,line=2.5,cex=1.2)

abline(h=sum(yu.p<0.05, na.rm=T)/sum(!is.na(yu.p)), col="gray")


dev.off()


