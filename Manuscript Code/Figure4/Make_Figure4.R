setwd("...")
source(".../LCO_generator_all_n.R")

load("test1_alt_all.Rdata")


# Load data for null 
null.table <- matrix(NA, nrow=1, ncol=4)
null.table[1,1] <- sum(yu.p<0.05)/length(yu.p)

load("test2_null_65new.Rdata")
load("test2_null_9new.Rdata")
load("test2_null_115new.Rdata")

null.table[1,2] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
null.table[1,3] <- sum(quantile.p.9.subs<0.05)/length(quantile.p.9.subs)
null.table[1,4] <- sum(quantile.p.115.subs<0.05)/length(quantile.p.115.subs)

# Run the LCO generator for a sample size of reps (1000 currently)
lco.bounds <- LCO.CI(reps,0.95,3)

bp.null <- boxplot(null.table)
bp.null$stats <- sapply(null.table,function(x) c(lco.bounds[(x*reps),2],x,x,x,lco.bounds[(x*reps),3]))


pdf(".../Fig4_null2.pdf",
    height=4.2,width=6.7)

plot.new()
plot.window(xlim=c(0.5,4.5),ylim=c(0,0.6),xaxs="i",yaxs="i")
rect(1.5,0,4.5,1,border="gray95",col="gray95")
bxp(bp.null,outline=F,medpch=16,medcex=0,whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.08,axes=F,add=T)
abline(h=0.05,col="gray60",lty=2)
axis(2,at=seq(0,0.6,by=0.1),label=c(0,0.1,0.2,0.3,0.4,0.5,0.6))
mtext("Estimated Type I Error rate",line=2.5,side=2,at=0.3)
axis(1,at=seq(1,4),label=c("YT06", "6.5A", "9A","11.5A"))
mtext("Null 2: SASA and compensatory processes", cex=1.5)
box()

dev.off()

