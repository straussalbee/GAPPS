
########## SET WD AND READ IN DATA #############################
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(scales)

load("freqTable.Rda") # frequencies from Flowjo
load("countsTable.Rda")#counts from Flowjo

#combine all data for renaming and organization
allData<-rbind(countsTable,freqTable)



######################## ACQUISITION NOTES #######################
#renaming the sample names: vaginal cells are all from donor 289. 
#"289" is the sample from 30Apr, "289b" is from 1May
#We had to switch out the nebulizer in the middle of the pbmc sample on
#1May so the first half of that sample is PBMC C and the 2nd half is
#PBMC D. I am treating them as different samples since we
#had to switch to a new nebulizer and couldn't 
#couldn't confirm the correct tuning for it since something was wrong
#with the tuning template.

#NOTE samples that are in dated fcs folders but not here were 
#removed from analysis because of low acquisition or other collection problems
#or because we didn't have anything to compare them to (i.e. facs tubes
#vs eppendorf pairings were incomplete)


########################## FIXING SAMPLE NAMES #########################
names(allData)[1]<-"SampleName"
allData<-allData%>%
  filter(SampleName!="Mean"&SampleName!="SD")%>%
  select(-X.1)
#add a column differentiating between count and freq data
allData<-mutate(allData, Stat = rep(c("Count","Freq"),times=c(11,11)))


#get rid of the brackets, escaping the bracket symbol so R doesn
#think that it is a regex group
allData$SampleName<-str_replace_all(allData$SampleName,"\\[ ","")
allData$SampleName<-str_replace_all(allData$SampleName," \\]","")



#################### ADD DATE COLLECTED ############################

#vector of dates, repeated times = # of samples on that date
Date <- rep(c("Apr29","Apr30","May1"),times = c(3,3,5))


#make a df with samples in order of date                       
DateCollected<-data.frame("SampleName"=c("B2_5mlEppendorf.fcs","hladik pbmc.fcs",
         "cytobrush 5084a-concat.fcs",
         "cytobrush 4897.fcs","hladik pbmc eppendorf.fcs",
         "vaginal 289.fcs","4884b.fcs","5168.fcs",
         "pbmc.fcs","pbmc2.fcs","vaginal.fcs"),
         "NewSampleName"=c("Blish PBMC","PBMC A","5084","4087","PBMC B",
                           "Vaginal A","4884 B","5168","PBMC C","PBMC D",
                           "Vaginal B"))

DateCollected$Date<-Date
#Merge the two data frames by Sample
allData<-merge(allData,DateCollected)


################## MAKE "OTHER" COLUMN FOR FREQ #############
#define a df with just frequency data and make a column for "Other" cell
#which = 100 - the sum of the frequency of the specified cell types
#note that CD4 and Cd8 frequencies are already included in Cd3 so they
#should be excluded for this calcuation. Also exclude AllLive because
#that is a parent of all gates below that
freq<-allData%>%
  filter(allData$Stat=="Freq")%>%
  mutate(Other = 100-Neutrophils-T.cells-NK-CD14)
               
freqMelt<-melt(freq, id.vars =c("NewSampleName","SampleName","Date",
                                "Stat"),
               variable.name="Cell.Type", value.name="Frequency")


freqMelt<-mutate(freqMelt, type = ifelse(str_detect(freqMelt$NewSampleName,
                                                    "PBMC")==TRUE,
                                         "PBMC","Cytobrush"))

freqMelt$type<-ifelse(str_detect(freqMelt$NewSampleName,"V")==TRUE,
                 "Vaginal",freqMelt$type)


######################### FREQUENCY PLOTS #########################
#Remove AllLive since it is not useful
freqMelt<-filter(freqMelt,Cell.Type!="AllLive")


#stacked bar
#remove Cd4 and cd8 so total = 100% (they are subsets of T cells)
freqMeltStackedBar<-filter(freqMelt, 
                            Cell.Type!= "CD4" & Cell.Type!="CD8")

#Leave out Blish PBMC since it was a bad sample

freqMeltStackedBar<-filter(freqMeltStackedBar,NewSampleName!="Blish PBMC")

#now plot
stackedBar<-ggplot(freqMeltStackedBar, aes(x=NewSampleName, y = Frequency))+
  geom_bar(aes(fill=Cell.Type),stat="identity")+
  labs(x="Sample")+
  scale_fill_manual(name="Cell Type",values=c("#1b9e77","#4daf4a","#984ea3",
                             "#e7298a","#80cdc1"))+
  theme(axis.text.x=element_text(angle=20,vjust=.8, size=10),
        axis.title=element_text(size= 15),
        legend.title=element_text(size=15),
        axis.title.x=element_text(vjust=.5),
        plot.title = element_text(vjust=2))+
  ggtitle("Frequency of (live) cell types in cytobrushes, PBMC and vaginal
cell samples")


#SAVE
ggsave(filename = "stackedBar.png", dpi=600)


#boxplot with points behind
#to remove points that are duplicated (outliers from boxplots)
#set outlier.color=NA

ggplot(freqMelt, aes(x=Cell.Type, y= Frequency))+
  geom_boxplot(outlier.colour=NA,fill=NA)+
  geom_point(aes(color=type),size=4)+
  labs(x="Cell Type",size=15)+
  ggtitle("Frequency of Cell Types in Cytobrushes, Vaginal cells \n\
and PBMC samples")+
  theme(axis.text.x=element_text(size=12),
        axis.title=element_text(size= 15))
#points
ggplot(freqMelt, aes(x=Cell.Type, y= Frequency))+
  geom_point(aes(color=type),size=5,alpha=0.8)

############################ COUNT DATA #############################
## create  a df for count data the same way as for freq data above
##remove Blish PBMC
count<-allData%>%
  filter(allData$Stat=="Count", NewSampleName!="Blish PBMC")

countMelt<-melt(count, id.vars =c("NewSampleName","SampleName","Date",
                                "Stat"),
               variable.name="Cell.Type", value.name="Count")


countMelt<-mutate(countMelt, type = ifelse(str_detect(countMelt$NewSampleName,
                                                    "PBMC")==TRUE,
                                         "PBMC","Cytobrush"))

countMelt$type<-ifelse(str_detect(countMelt$NewSampleName,"V")==TRUE,
                      "Vaginal",countMelt$type)


########################## COUNT TABLES #################################
library(pander)
NeutrophilCounts<-countMelt%>%
  filter(Cell.Type %in% c("Neutrophils","AllLive"))%>%
  arrange(type,Cell.Type)%>%
  select(NewSampleName,Count,Cell.Type,type)




NeutrophilCountsWide<-dcast(NeutrophilCounts, NewSampleName~Cell.Type,
      value.var = "Count")

pander(NeutrophilCountsWide)


#Change the sample names a little to make it look better on the plot

NeutrophilCounts$NewSampleName<-str_replace_all(NeutrophilCounts$NewSampleName,
                                                "PBMC","")
###################### COUNT PLOTS ################################

NeutroAllCount<-ggplot(NeutrophilCounts,aes(x = NewSampleName, y = Count))+
  geom_point(aes(color=Cell.Type),size=4, alpha=0.6)+
  labs(x = "Sample")+
  scale_y_continuous(labels=comma)+
  scale_color_discrete(name="Cells")+
  ggtitle("Count of Neutrophils and total live cells in cytobrushes,PBMC and
vaginal cell samples")+
  facet_wrap(~type, scales="free")+
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size= 15),
        axis.title.x=element_text(vjust=.5),
        legend.title=element_text(size=15),
        plot.title = element_text(vjust=2),
        strip.text.x=element_text(size=14),
        aspect.ratio=1.5)

ggsave("NeutroAllCount.png",dpi=600)
