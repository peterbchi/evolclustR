
setwd("...")
sens <- NA
mspc <- NA

load("test2_surface2_rads.RData")
sens[1] <- sum(yu.p<0.05)/length(yu.p)
sens[4] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
load("test2_surface3_rads.Rdata")
sens[2] <- sum(yu.p<0.05)/length(yu.p)
sens[5] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
load("test2_surface4_rads.Rdata")
sens[3] <- sum(yu.p<0.05)/length(yu.p)
sens[6] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)

load("test1_alt_all.Rdata")
mspc[1] <- sum(yu.p<0.05)/length(yu.p)
mspc[2] <- sum(yu.p<0.05)/length(yu.p)
mspc[3] <- sum(yu.p<0.05)/length(yu.p)
load("test2_null_65new.Rdata")
mspc[4] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
mspc[5] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
mspc[6] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs) - 0.005


pdf(".../Fig5_roc.pdf",
    height=6,width=6)

plot(sens~mspc,xlim=c(0,1), ylim=c(0,1),pch=c(17,17,17,19,19,19),
     col=c("gray10", "gray35","gray70"),
     xlab="Type I Error rate", ylab="Power",main="")
legend("bottomright",c("YT06, 2 subs", "YT06, 3 subs", "YT06, 4 subs",
                       "CKL20, 2 subs", "CKL20, 3 subs", "CKL20, 4 subs"),
       pch=c(17,17,17,19,19,19), col=c("gray10", "gray35","gray70"))
abline(a=0,b=1,lty=2)
dev.off()
