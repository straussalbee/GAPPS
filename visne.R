library(flowCore)
library(Rtsne)
library(ggplot2)
library(reshape2)
MICC_15_130_01 <- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\15-130-01_normalized.exported.FCS3.fcs")
HVTN5217 <- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\5217.exported.FCS3.fcs")
HVTN4396<- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\4396_normalized.exported.FCS3.fcs")
HVTN3633<-read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\3633_normalized.exported.FCS3.fcs")

cellsPerSample <- 300

# get random sample of cells from each FCS file
MICC_15_130_01_sub <- MICC_15_130_01[sample(1:nrow(MICC_15_130_01), cellsPerSample), ]
HVTN5217_sub <- HVTN5217[sample(1:nrow(HVTN5217), cellsPerSample), ]
HVTN4396_sub <- HVTN4396[sample(1:nrow(HVTN4396), cellsPerSample), ]
HVTN3633_sub <- HVTN3633[sample(1:nrow(HVTN3633), cellsPerSample), ]

# combine files
# skip columns 1 and 2 (time and event length)(S4)
# @exprs removes matrix from flowframe
#one row per cell, columns are channels
#rbind multiple samples
toCluster <- rbind(MICC_15_130_01_sub@exprs[,3:50],HVTN5217_sub@exprs[,3:50],
                   HVTN4396_sub@exprs[,3:50],HVTN3633_sub@exprs[ , 3:50])#leaves out cell length and time

# arcsinh transformation (something about cofactor of 5 in visne paper that I don't understand)
toCluster <- asinh(toCluster)

# Generating viSNE maps included the following steps (exact details can be 
# found in Supplementary Table 1). First, between 6,000 and 12,000 cells were 
# uniformly subsampled from the data. After subsampling, viSNE was run for 500 
# iterations to project the data into 2D. Unlike t-SNE, PCA was not used as a 
# preprocessing step. All runs used an identical random seed and the default 
# t-SNE parameters (perplexity = 30, momentum = 0.5 for initial 250 iterations, 
# momentum = 0.8 for remaining iterations, epsilon = 500, lie factor = 4 for 
# initial 100 iterations, lie factor = 1 for remaining iterations).
set.seed(123) # always get same result with same data
cluster <- Rtsne(toCluster, initial_dims = 10^6, max_iter = 500, PCA = FALSE, 
   perplexity = 30, theta = 0)

# tsne1 and tsne2 are in the matrix in cluster$Y
colnames(cluster$Y) <- c("tsne1", "tsne2")

# combine tsne dimensions with the mass cytometry data
toPlot <- cbind(toCluster, cluster$Y)
toPlot <- as.data.frame(toPlot)

# label samples
toPlot$Sample <- c(rep("MICC_15_130_01", cellsPerSample),
                   rep("HVTN5217", cellsPerSample),
                   rep("HVTN4396", cellsPerSample),
                   rep("HVTN3633", cellsPerSample))



myPlot<-ggplot(as.data.frame(toPlot), aes(tsne1, tsne2)) +
   geom_point(aes(color = Sample)) 
#+scale_color_gradientn(colours=topo.colors(10))
ggsave("myPlot.png",dpi=600)

meltToPlot<-melt(toPlot,id.vars="Sample")
names(meltToPlot)[3]<-"channel"

ggplot(meltToPlot, aes(channel),position="dodge") +
  geom_bar(aes(fill = variable))
