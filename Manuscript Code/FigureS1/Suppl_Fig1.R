setwd("...")
load(".../test2_surface2_rads.RData")

alt.table <- matrix(NA, nrow=1, ncol=7)
alt.table[1,1] <- sum(yu.p<0.05)/length(yu.p)
alt.table[1,2] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
alt.table[1,3] <- sum(quantile.p.7.subs<0.05)/length(quantile.p.65.subs)
alt.table[1,4] <- sum(quantile.p.75.subs<0.05)/length(quantile.p.65.subs)
alt.table[1,5] <- sum(quantile.p.8.subs<0.05)/length(quantile.p.65.subs)
alt.table[1,6] <- sum(quantile.p.85.subs<0.05)/length(quantile.p.65.subs)
alt.table[1,7] <- sum(quantile.p.9.subs<0.05)/length(quantile.p.9.subs)



load(".../test2_surface2_rads_1d4t.RData")
alt.table1d4t <- matrix(NA, nrow=1, ncol=13)
alt.table1d4t[1,1] <- sum(yu.p<0.05)/length(yu.p)
alt.table1d4t[1,2] <- sum(quantile.p.45.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,3] <- sum(quantile.p.5.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,4] <- sum(quantile.p.55.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,5] <- sum(quantile.p.6.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,6] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,7] <- sum(quantile.p.7.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,8] <- sum(quantile.p.75.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,9] <- sum(quantile.p.8.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,10] <- sum(quantile.p.85.subs<0.05)/length(quantile.p.65.subs)
alt.table1d4t[1,11] <- sum(quantile.p.9.subs<0.05)/length(quantile.p.9.subs)
alt.table1d4t[1,12] <- sum(quantile.p.95.subs<0.05)/length(quantile.p.95.subs)
alt.table1d4t[1,13] <- sum(quantile.p.10.subs<0.05)/length(quantile.p.10.subs)

load(".../test2_surface2_rads_1AX8.RData")
alt.table1ax8 <- matrix(NA, nrow=1, ncol=13)
alt.table1ax8[1,1] <- sum(yu.p<0.05)/length(yu.p)
alt.table1ax8[1,2] <- sum(quantile.p.45.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,3] <- sum(quantile.p.5.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,4] <- sum(quantile.p.55.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,5] <- sum(quantile.p.6.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,6] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,7] <- sum(quantile.p.7.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,8] <- sum(quantile.p.75.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,9] <- sum(quantile.p.8.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,10] <- sum(quantile.p.85.subs<0.05)/length(quantile.p.65.subs)
alt.table1ax8[1,11] <- sum(quantile.p.9.subs<0.05)/length(quantile.p.9.subs)
alt.table1ax8[1,12] <- sum(quantile.p.95.subs<0.05)/length(quantile.p.95.subs)
alt.table1ax8[1,13] <- sum(quantile.p.10.subs<0.05)/length(quantile.p.10.subs)


pdf(".../Suppl_Fig1.pdf")
par(mfrow=c(3,1),oma=c(3,2.5,1.5,1),mar=c(3,2,1.5,1))

radii <- seq(6.5,9,by=0.5)
plot(alt.table[2:7]~radii, type='l',ylim=c(0,1),ylab="",xlab="",xlim=c(6.5,9),main="2I0Q")
abline(h=alt.table[1], col="gray")
mtext("power",side=2,line=2.5,cex=1.2)
legend("topright", legend=c("CKL20","YT06"),col=c("black","gray"), lty=1, cex=1.5)


plot(alt.table1d4t[6:11]~radii, type='l',ylim=c(0,1),ylab="Power",xlab="",main="1D4T")
abline(h=alt.table1d4t[1], col="gray")
mtext("power",side=2,line=2.5,cex=1.2)



plot(alt.table1ax8[6:11]~radii, type='l',ylim=c(0,1),ylab="Power", main="1AX8")
abline(h=alt.table1ax8[1], col="gray")
mtext("power",side=2,line=2.5,cex=1.2)
mtext("radius",side=1,line=3,cex=1.2)

dev.off()


