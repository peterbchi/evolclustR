setwd("...")
source(".../LCO_generator_all_n.R")


load("test1_null_all.RData")


# Load data for null 
null.table <- matrix(NA, nrow=1, ncol=4)
null.table[1,1] <- sum(yu.p<0.05)/length(yu.p) # need to run sims to get this
null.table[1,2] <- sum(quantile.p.65.subs<0.05)/length(quantile.p.65.subs)
null.table[1,3] <- sum(quantile.p.9.subs<0.05)/length(quantile.p.9.subs)
null.table[1,4] <- sum(quantile.p.115.subs<0.05)/length(quantile.p.115.subs)

# Run the LCO generator for a sample size of reps (1000 currently)
lco.bounds <- LCO.CI(reps,0.95,3)

bp.null <- boxplot(null.table)
bp.null$stats <- sapply(null.table,function(x) c(lco.bounds[(x*reps),2],x,x,x,lco.bounds[(x*reps),3]))


pdf(".../Fig2_null1.pdf",
    height=4.2,width=6.7)

plot.new()
plot.window(xlim=c(0.5,4.5),ylim=c(0,0.15),xaxs="i",yaxs="i")
rect(1.5,0,4.5,1,border="gray95",col="gray95")
bxp(bp.null,outline=F,medpch=16,medcex=0,whisklty=1,whisklwd=1.5,staplewex=1,staplelwd=1.5,boxwex=0.08,axes=F,add=T)
abline(h=0.05,col="gray60",lty=2)
axis(2,at=seq(0,0.15,by=0.05),label=c(0,0.05,0.1,0.15))
mtext("Estimated Type I Error rate",line=2.5,side=2,at=0.075)
axis(1,at=seq(1,4),label=c("YT06", "6.5A", "9A","11.5A"))
mtext("Simulations under naive null hypothesis", cex=1.5)
box()

dev.off()

