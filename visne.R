library(flowCore)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)

################## READ IN FCS FILES (LIVE POPN, EXPORTED) ###########

#locate the folder where the fcs files are 
folderPath<-"J:\\MacLabUsers\\Claire\\Projects\\GAPPS Project\\GAPPS 2015 Neutrophils\\data and analysis\\fcs and compiled flowjo\\normed live export"

#use this later to name the flow frames
files<-list.files(folderPath)
fileNamesOnly<-as.list(str_replace(files,pattern=".fcs",replacement=""))


#The flowCore method for reading in files
frames<-lapply(dir(folderPath, pattern = "\\.fcs",full.names = TRUE),
               read.FCS)


as(frames,"flowSet") #check

#name the list of frames with the lists of file names(minus the .fcs)
names(frames) <- fileNamesOnly

fs<-as(frames,"flowSet") #now I have a flowSet of multiple flowFrames


#extract the data matrix from the flowSet to check if it looks right
matrix<-fsApply(fs,exprs)

#the first and last 6 rows of this is the same matrix as
#if you read in the files just as a LIST from the folder 
#(lappy(files,read.FCS)) and then do lapply(data,exprs)
# instead of using the flowSet methods so I think both methods
#get you the same result, i.e. a single data matrix that includes
#data from all the samples together, in the order that they
# were read into the flowSet.

#Note: if one wants to compare the matrix you get from reading in
#the files in different ways, you need to do it before the random
#subsampling. But maybe if you do set.seed before the random sampling
#it will end up the same? I haven't tried this yet 24Aug15




###### SUBSAMPLE THE SAMPLES #####################

cellsPerSample <- 500

# here is a function to get random sample of cells from each FCS file

subSample<-function(Frame){
  Frame[sample(1:nrow(Frame),cellsPerSample),]
}

#run the function on the flowSet
sampledFs<-fsApply(fs,FUN=subSample)




########################### COMBINE FILES ###########################

# arcsinh transformation with cofactor of 5 (divide data by 5)

#here is a function to get the actual data matrix out of the subsampled
#flowframe
# exprs() removes matrix from flowframe
#one row per cell, columns are channels
getMatrix<-function(subSampledSet){
  exprs(subSampledSet)
}

# apply the function to the subSampled flowset
 toCluster<-fsApply(sampledFs,FUN=getMatrix)

toCluster<-asinh(toCluster/5)# this is what the nolan lab does
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

#This next step takes a minute to run...

cluster<-Rtsne(toCluster, initial_dims = 10^6, max_iter = 500, PCA = FALSE, 
                perplexity = 30, theta = 0)


# tsne1 and tsne2 values are in the matrix in cluster$Y
#rename them for clarity

colnames(cluster$Y) <- c("tsne1", "tsne2")


#combine tsne dimensions from the Rtsne function
#with the mass cytometry data extracted and transformed from the .fcs files

toPlot <- cbind(toCluster, cluster$Y)
toPlot <- as.data.frame(toPlot)

#toPlot is a df with one column per metal and one row per sample

# label samples: make a column in the toPlot df called
#Sample and repeat the sample name for the number cells sampled

toPlot$Sample <- c(rep("MICC_15_130_01", cellsPerSample),
                   rep("HVTN5217", cellsPerSample),
                   rep("HVTN4396", cellsPerSample),
                   rep("HVTN3633", cellsPerSample))
head(toPlot)

#shows each sample as a different color
tSNE4samples<-ggplot(as.data.frame(toPlot), aes(tsne1, tsne2)) +
  geom_point(aes(color = Sample)) 
ggsave("tSNE4samples10Aug15.png",dpi=600)


######################  GET CHANNEL NAMES and MELT  ########################
#I am doing this so I can make gradient plots of tsne1 by tsne2 for
#each Ab and have the plot labeled with the Ab instead of with the
#name of the metal is is conjugated to.

#extract the names assigned to the channels from the flowframe(ex. CD3)
#they are the same for all samples so I can get them from any of the 
#flow frames

pData<-as.data.frame(pData(parameters(frames[[1]])))

#select just the "name" (metal) and the "desc" (the ab assigned to it)
channelAndDesc<-select(pData,name,desc)

#changing to characters instead of factors
channelAndDesc$name<-as.character(channelAndDesc$name)

#fix the places where there is no entry for desc
#i.e. channel name is Event_length and Time

channelAndDesc$desc[channelAndDesc$name=="Event_length"]<-"Event Length"
channelAndDesc$desc[channelAndDesc$name=="Time"]<-"Time"

channelAndDesc$desc<-as.character(channelAndDesc$desc)

#get rid of the rownames (which indicate keywords??)
rownames(channelAndDesc)<-NULL

######## MELT #####################################
#now melt the df so you can plot tsne1 vs tsne2 (like in viSNE)and

meltToPlot<-melt(toPlot,id.vars=c("tsne1","tsne2","Sample"))

###################### MERGE IN NAMES #######################
#merge in the df of the channels and their corresponding names
#so they Ab names can be used in the plot

meltToPlot<-merge(meltToPlot, channelAndDesc,by.x="variable",by.y="name")

allParams<-ggplot(meltToPlot, aes(x=tsne1,y=tsne2))+
  geom_point(aes(color = value),alpha=0.6, size=.6) +
  scale_colour_gradientn("value",colours=topo.colors(7))+
  facet_wrap(~desc) +
  xlab("bhSNE1") + 
  ylab("bhSNE2") +
    #scale_x_continuous(limits=c(-40,40)) +
    #scale_y_continuous(limits=c(-40,40)) +
  theme(panel.grid.major=element_blank())

ggsave("allParams10Aug15.png",dpi=600)
                           
#for scale colour gradientn("value, colours=topo.colors(7))