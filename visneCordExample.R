#150519 plot viSNE cord data 
require(ggplot2)
require(reshape2)
require(plyr)
require(gridExtra)
require(stringr)

#final version with emily
setwd("J:/MacLabUsers/Claire/Projects/GAPPS Project/GAPPS 2015 Neutrophils/data")

# Read in and subset viSNE data
#cordvisne <- read.csv("150519 viSNE Cord.csv", header=TRUE)
cordvisne <- read.csv("f.csv", header=TRUE)
head(cordvisne)

names(cordvisne)<-gsub('.*_',"",names(cordvisne))

names(cordvisne)<-substr(names(cordvisne),1,nchar(names(cordvisne))-1)
#mostly right

cordvisne.subset <- subset(cordvisne, 
                           select=c("CD45","CD66b","CD4","CD8","CD3","CD16","Cd15","X.bh.SNE1",
                                    "X.bh.SNE2"))
head(cordvisne.subset)
cordvisne.melt <- melt(cordvisne.subset,id=c("X.bh.SNE1", "X.bh.SNE2"))
names(cordvisne.melt)[3] <- "Marker"
head(cordvisne.melt)


adultvisne <- read.csv("viSNE Adult.csv", header=TRUE)
head(adultvisne)
adultvisne.subset <- subset(adultvisne, 
                            select=c("bh.SNE1", "bh.SNE2", "CD57", "NKG2A", "Syk", "DNAM1", "CD8", 
                                     "NKp30", "NKp46", "CD16", "K3DL1", "K2DL1", "K2DL3", "CD56", 
                                     "PD1", "A2B4", "K2DS4","LILRB1", "FcRg", "NKG2C"))
head(adultvisne.subset)
adultvisne.melt <- melt(adultvisne.subset,id=c("bh.SNE1", "bh.SNE2"))
names(adultvisne.melt)[3] <- "Marker"
head(adultvisne.melt)

# Plot data
cordPlots <- ggplot(cordvisne.melt, aes(X.bh.SNE1, X.bh.SNE2)) + 
  geom_density2d(colour="gray50", size = 0.2) +
  geom_point(aes(colour=value), alpha=0.6, size=0.5) +
  scale_colour_gradientn("Value", colours=topo.colors(7), limits=c(-1,6.5)) +
  facet_wrap(~ Marker) +
  xlab("bhSNE1") +
  ylab("bhSNE2") +
  #scale_x_continuous(limits=c(-40,40)) +
  #scale_y_continuous(limits=c(-40,40)) +
  theme(panel.grid.major=element_blank())
cordPlots
ggsave("cordvisne.pdf", width=12, height=10)


adultPlots <- ggplot(adultvisne.melt, aes(bh.SNE1, bh.SNE2)) + 
  geom_density2d(colour="gray50", size = 0.2) +
  geom_point(aes(colour=value), alpha=0.6, size=0.5) +
  scale_colour_gradientn("Value", colours=topo.colors(7), limits=c(-1,6.5)) +
  facet_wrap(~ Marker) + 
  xlab("bhSNE1") +
  ylab("bhSNE2") +
  scale_x_continuous(limits=c(-40,40)) +
  scale_y_continuous(limits=c(-40,40)) +
  theme(panel.grid.major=element_blank())
adultPlots
ggsave("adultvisne.pdf", width=12, height=10)