setwd("C:/Users/pchi01/Dropbox/Temple/Clustering/MBE-Resubmission/leptin/Figure6")
library(tidyverse)
library(gridExtra)

load("test2_surface3_hominoid_1ax8.RData")

df <- data.frame(
  data_type = factor(c(rep("Simulated", length(as.vector(dists_1ax8))), rep("Real Data", length(dists_hominoid_1ax8)))),
                 dist = c(as.vector(dists_1ax8), dists_hominoid_1ax8)
)

plot1 <- ggplot(df, aes(x=dist, fill=data_type, y = ..density..)) + 
  geom_histogram(alpha=0.5, position="identity", bins=30) + 
  scale_fill_grey() +
  scale_x_continuous(name="pairwise distance") + 
  scale_y_continuous(limits=c(0,0.25)) +
  theme(legend.position = "none") + 
  ggtitle("Hominoid Branch")

#ggplot() + 
#  geom_histogram(aes(x =as.vector(dists_1ax8), y= ..density..), fill="gray20", alpha=0.2) +
#  geom_histogram(aes(x =as.vector(dists_hominoid_1ax8), y= ..density..), bins=15, fill="black", alpha=0.6) +
#  scale_x_continuous(name="pairwise distance", limits=c(0,50)) +
#  scale_color_manual(values=c("gray20", "black"))

load("test2_surface3_macaca_1ax8.RData")

df <- data.frame(
  data_type = factor(c(rep("Simulated", length(as.vector(dists_1ax8))), rep("Real Data", length(dists_hominoid_1ax8)))),
  dist = c(as.vector(dists_1ax8), dists_hominoid_1ax8)
)  # forgot to rename hominoid to macaca in script but it's right

plot2 <- ggplot(df, aes(x=dist, fill=data_type, y = ..density..)) + 
  geom_histogram(alpha=0.5, position="identity", bins=30) + 
  scale_fill_grey(name="data type") +
  scale_x_continuous(name="pairwise distance") +   
  scale_y_continuous(limits=c(0,0.25)) +
  ggtitle("Macaca Branch")

grid.arrange(plot1, plot2, nrow=1)
