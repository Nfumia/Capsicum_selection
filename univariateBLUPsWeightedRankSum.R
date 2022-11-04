### Weighted Rank Sum Index 
# derived from adjusting BLUPs by accuracy of the model predicting the BLUP, thereby reducing the weight of importance per phenotype

WorkDir <- "C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel"
setwd(WorkDir)

# read in .csv files of univariate model prediction accuracies and genoxpheno BLUPs
CorrGControl <- read.csv("univariateGBLUP/traitGCorrelations_KCV_control.csv",row.names=1)
BLUPGControl <- read.csv("univariateGBLUP/traitGBLUPS_KCV_control.csv",row.names=1)

CorrGStress1 <- read.csv("univariateGBLUP/traitGCorrelations_KCV_stress1.csv",row.names=1)
BLUPGStress1 <- read.csv("univariateGBLUP/traitGBLUPS_KCV_stress1.csv",row.names=1)

CorrGStress2 <- read.csv("univariateGBLUP/traitGCorrelations_KCV_stress2.csv",row.names=1)
BLUPGStress2 <- read.csv("univariateGBLUP/traitGBLUPS_KCV_stress2.csv",row.names=1)

phenos <- read.csv("univariatePHENOBLUP/traitBLUPs_combined.csv")

# data -> numeric
BLUPGControl[1:75] <- sapply(BLUPGControl[1:75],as.numeric)
BLUPGStress1[1:75] <- sapply(BLUPGStress1[1:75],as.numeric)
BLUPGStress2[1:75] <- sapply(BLUPGStress2[1:75],as.numeric)

# normalize data between [0,1] (equal scale of phenotypes for index)
BLUPGControlscale <- sapply(BLUPGControl[1:75], function(x) (x-min(x))/(max(x)-min(x)))
rownames(BLUPGControlscale) <- rownames(BLUPGControl)
BLUPGControlscale <- as.data.frame(BLUPGControlscale)
BLUPGStress1scale <- sapply(BLUPGStress1[1:75], function(x) (x-min(x))/(max(x)-min(x)))
rownames(BLUPGStress1scale) <- rownames(BLUPGStress1)
BLUPGStress1scale <- as.data.frame(BLUPGStress1scale)
BLUPGStress2scale <- sapply(BLUPGStress2[1:75], function(x) (x-min(x))/(max(x)-min(x)))
rownames(BLUPGStress2scale) <- rownames(BLUPGStress2)
BLUPGStress2scale <- as.data.frame(BLUPGStress2scale)

# extract ranking of genotypes per phenotype

library(dplyr)
datC <- BLUPGControlscale %>%
  mutate(across(c(1:75), rank))

datS1 <- BLUPGStress1scale %>%
  mutate(across(c(1:75), rank))

datS2 <- BLUPGStress2scale %>%
  mutate(across(c(1:75), rank))


# weight rankings per phenotype by the univariate accuracy of prediction on that trait 
# ranking * (1-model accuracy) => inflates importance of a trait (by reducing rank which is good) proportionally to confidence in model

for (i in 1:ncol(datC)) {
  num <- CorrGControl[i,]
  datC[,i] <- datC[,i]*(1-num)
}

for (i in 1:ncol(datS1)) {
  num <- CorrGStress1[i,]
  datS1[,i] <- datS1[,i]*(1-num)
}

for (i in 1:ncol(datS2)) {
  num <- CorrGStress2[i,]
  datS2[,i] <- datS2[,i]*(1-num)
}


# sum rankings for phenotypes of each genotypes (lowest number is best!)
datC <- datC %>% 
  replace(is.na(.), 0) %>%
  mutate(weightRankSum=rowSums(.[,-1]))

datS1 <- datS1 %>% 
  replace(is.na(.), 0) %>%
  mutate(weightRankSum=rowSums(.[,-1]))

datS2 <- datS2 %>% 
  replace(is.na(.), 0) %>%
  mutate(weightRankSum=rowSums(.[,-1]))

# pull genotype name and rank sum

genoRankSumC <- datC[76]
genoRankSumC <- genoRankSumC[with(genoRankSumC,order(weightRankSum)),]
genoRankSumC <- as.data.frame(genoRankSumC)
genoRankSumC$g2pname <- rownames(genoRankSumC)

genoRankSumS1 <- datS1[76]
genoRankSumS1 <- genoRankSumS1[with(genoRankSumS1,order(weightRankSum)),]
genoRankSumS1 <- as.data.frame(genoRankSumS1)
genoRankSumS1$g2pname <- rownames(genoRankSumS1)

genoRankSumS2 <- datS2[76]
genoRankSumS2 <- genoRankSumS2[with(genoRankSumS2,order(weightRankSum)),]
genoRankSumS2 <- as.data.frame(genoRankSumS2)
genoRankSumS2$g2pname <- rownames(genoRankSumS2)

# merge the 3 Rank Sums (1 for each environment)
mergeRankSum <- merge(genoRankSumC,genoRankSumS1,by="g2pname",all=T) %>% 
    merge(genoRankSumS2)
str(mergeRankSum)

mergeRankSum <- mergeRankSum %>%
  mutate(across(c(2:4), rank))

mergeRankSum$genotype <- phenos$genotype[match(mergeRankSum$g2pname,phenos$g2pname)]

mergeRankSum$meanRank <- rowMeans(mergeRankSum[2:4])

mergeRankSum$meanRankStress <- rowMeans(mergeRankSum[3:4])

#write.csv(mergeRankSum,"C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel/WeightedRankSum.csv",row.names=FALSE)

# column of rank sum deviation with increasing temperature environments 
mergeRankSum$deviationRankSum <- abs(mergeRankSum$genoRankSumC-mergeRankSum$genoRankSumS1-mergeRankSum$genoRankSumS2)


df <- subset(mergeRankSum[,c(1:7)])

phenosMean <- phenos %>% 
  group_by(g2pname) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

df <- merge(df,phenosMean,by.x="g2pname")
df[,c(8:82)] <- scale(df[,c(8:82)])

library(tidyr)
long <- df %>% gather(key=rankID,value=rank,genoRankSumC:genoRankSumS2,-c(1,5:82))
long <- long %>% gather(key=trait,value=BLUP,DASanthesis_BLUP:biomassAUCperplant_BLUP,-c(1,2))
long$trait<-gsub("_BLUP","",as.character(long$trait))

long <- subset(long,trait!="biomassAUC")
long <- subset(long,trait!="biomassAUCperplant")

long$rankID[long$rankID=="genoRankSumC"] <- "Control"
long$rankID[long$rankID=="genoRankSumS1"] <- "Increased Temp"
long$rankID[long$rankID=="genoRankSumS2"] <- "Severe Temp"


#long$logBLUP <- log(long$BLUP)

library(ggplot2)
library(ggforce)
p<- ggplot(data=long,aes(x=BLUP,y=rank,group=trait,color=rankID)) +
  geom_point(aes(x=BLUP,y=rank,color=rankID),size=.02) +
  geom_smooth(data=long,aes(group=rankID,color=rankID),method="lm",se=FALSE,formula=y~x,size=.5) +
  facet_wrap_paginate(facets=~trait,scales="free",nrow=4,ncol=4,page=5) +
  labs(color="Environment")
ggsave('ranksumXtraits_plot5.pdf',plot=p,width=8.5,height=11,device="pdf",units="in")




  


with(df,plot(ndvimean_BLUP,genoRankSumC,xlab="Mean NDVI",ylab="Weighted Rank Sum Index",main="Correlation of Index and Trait by Environment"))
abline(lm(df$genoRankSumC~df$ndvimean_BLUP))
points(df$ndvimean_BLUP,df$genoRankSumS1,col="red")
abline(lm(df$genoRankSumS1~df$ndvimean_BLUP),col="red")
points(df$ndvimean_BLUP,df$genoRankSumS2,col="blue")
abline(lm(df$genoRankSumS2~df$ndvimean_BLUP),col="blue")
legend(0.705,275,legend=c("Control","High Temp","Severe Temp"),
       col=c("black","red","blue"), lty=1, cex=0.8)
