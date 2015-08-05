library(flowCore)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(dplyr)

################## READ IN FCS FILES (LIVE POPN, EXPORTED) ###########



MICC_15_130_01 <- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\15-130-01_normalized.exported.FCS3.fcs")
#HVTN5217 <- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\5217.exported.FCS3.fcs")
HVTN4396<- read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\4396_normalized.exported.FCS3.fcs")
HVTN3633<-read.FCS("J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\for viSNE\\3633_normalized.exported.FCS3.fcs")

#checking the parameters
pDataMICC<-pData(parameters(MICC_15_130_01))
pData5217<-pData(parameters(HVTN5217))
pData4396<-pData(parameters(HVTN4396))
pData3633<-pData(parameters(HVTN3633))

identical(pData3633$name,pData4396$name)
identical(pData3633$name,pDataMICC$name)

identical(pData3633$desc,pDataMICC$desc)
identical(pData3633$desc,pData4396$desc)

#which names that are in 5217 are NOT in 3633
subset(pData5217$name,!(pData5217$name %in% pData3633$name))
#$P49N     $P50N 
#"Cd106Di" "Cd108Di" 
#these are extra qdot channels, also the desc and name aren't matching up correctly

################# ISSUE ######################################
#5Aug15: Found issue when trying to remove unwanted channels from
#the FCS files that I read in, some names and desc don't match up.
#I checked flowjo and if I open the exported file in a new wsp I see
#the same problem. Things match ok in the original wsp.
#Tried to re-export the live population 2x and got "export error"
#leaving that sample out for now


###### SUBSAMPLE THE SAMPLES #####################

cellsPerSample <- 500

# get random sample of cells from each FCS file

MICC_15_130_01_sub <- MICC_15_130_01[sample(1:nrow(MICC_15_130_01), cellsPerSample), ]
#HVTN5217_sub <- HVTN5217[sample(1:nrow(HVTN5217), cellsPerSample), ]
HVTN4396_sub <- HVTN4396[sample(1:nrow(HVTN4396), cellsPerSample), ]
HVTN3633_sub <- HVTN3633[sample(1:nrow(HVTN3633), cellsPerSample), ]

######################  GET CHANNEL NAMES ########################
#extract the names assigned to the channels from the flowframe(CD3)
#they are the same for all samples so I can get them from any of the 
#flow frames

pData<-pData(parameters(HVTN3633_sub))

#select just the "name" (metal) and the "desc" (the ab assigned to it)
channelAndDesc<-select(pData,name,desc)
#changing to characters instead of factos
channelAndDesc$name<-as.character(channelAndDesc$name)
channelAndDesc$desc<-as.character(channelAndDesc$desc)

#get rid of the rownames (which indicate keywords??)
rownames(channelAndDesc)<-NULL


########################### COMBINE FILES ###########################
#I want to remove any columns that aren't relevant for visualization
#these are different for sample 5217 since I didn't collect beads
#For the other samples:
# skip columns 1,2,6,9,10,13,51,52 (time,event length,beads1, beads2,
# beadDist and Event #)
#Ce140 is beads1 amd Ce142 is beads2
#qdot4 = Cd113 qdot6=Cd116


# @exprs removes matrix from flowframe
#one row per cell, columns are channels
#rbind multiple samples
#Note that HVTN5217 doesn't have any of the beads columns so the indices to remove are different

toCluster <- rbind(MICC_15_130_01_sub@exprs[,-c(1,2,6,9,10,13,51,52)],
                   #HVTN5217_sub@exprs[,-c(1,2,6,9,51)],
                   HVTN4396_sub@exprs[,-c(1,2,6,9,10,13,51,52)],
                   HVTN3633_sub@exprs[ , -c(1,2,6,9,10,13,51,52)])
# arcsinh transformation with cofactr of 5 (divide data by 5)


toCluster2<-asinh(toCluster/5)# this is what the nolan lab does
#and probably viSNE because their default is cofactor of 5
#"xxxx2" going forward comes from the code above


# Generating viSNE maps included the following steps (exact details can be 
# found in Supplementary Table 1). First, between 6,000 and 12,000 cells were 
# uniformly subsampled from the data. After subsampling, viSNE was run for 500 
# iterations to project the data into 2D. Unlike t-SNE, PCA was not used as a 
# preprocessing step. All runs used an identical random seed and the default 
# t-SNE parameters (perplexity = 30, momentum = 0.5 for initial 250 iterations, 
# momentum = 0.8 for remaining iterations, epsilon = 500, lie factor = 4 for 
# initial 100 iterations, lie factor = 1 for remaining iterations).



set.seed(123) # always get same result with same data


#Rtsne gives a LIST of 7 attributes and their values
#theta, perplexity, N, origD,Y, costs and itercosts

cluster2<-Rtsne(toCluster2, initial_dims = 10^6, max_iter = 500, PCA = FALSE, 
                perplexity = 30, theta = 0)


# tsne1 and tsne2 values are in the matrix in cluster$Y

colnames(cluster2$Y) <- c("tsne1", "tsne2")


#combine tsne dimensions from the Rtsne functio
#with the mass cytometry data extracted and transformed from the .fcs files

toPlot2 <- cbind(toCluster2, cluster2$Y)
toPlot2 <- as.data.frame(toPlot2)


# label samples: make a column in the toPlot2 df called
#Sample and repeat the sample name for the number cells sampled

toPlot2$Sample <- c(rep("MICC_15_130_01", cellsPerSample),
                   #rep("HVTN5217", cellsPerSample),
                   rep("HVTN4396", cellsPerSample),
                   rep("HVTN3633", cellsPerSample))

#toPlot2 is a df with one column per metal and one row per sample

head(toPlot2)

#shows each sample as a different color
tSNE3samples<-ggplot(as.data.frame(toPlot2), aes(tsne1, tsne2)) +
  geom_point(aes(color = Sample)) 
ggsave("tSNE3samples.png",dpi=600)

#now melt the df so you can plot tsne1 vs tsne2 (like in viSNE)and
#color the plots by the value given for 
meltToPlot2<-melt(toPlot2,id.vars=c("tsne1","tsne2","Sample"))


#merge in the df of the channels and their corresponding names

meltToPlot2<-merge(meltToPlot2, channelAndDesc,by.x="variable",by.y="name")

allParams<-ggplot(meltToPlot2, aes(x=tsne1,y=tsne2))+
  geom_point(aes(color = value),alpha=0.6, size=.6) +
  scale_colour_gradientn("value", colours=topo.colors(7))+
  facet_wrap(~desc) +
  xlab("bhSNE1") + 
  ylab("bhSNE2") +
    #scale_x_continuous(limits=c(-40,40)) +
    #scale_y_continuous(limits=c(-40,40)) +
  theme(panel.grid.major=element_blank())

ggsave("allParams.png",dpi=600)
                           
