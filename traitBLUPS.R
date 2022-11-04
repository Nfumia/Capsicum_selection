### Load libraries, functions, and data ----
# Set directory where you have the data files as working dir.. 
WorkDir <- "C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel"
setwd(WorkDir)

# Import phenotypic data
pheno <- read.csv("combined.csv")
pheno[6:80] <- sapply(pheno[6:80],as.numeric)


### RUN LOOP WITH FULL DATA
phenoBLUP <- list()

for (i in 6:ncol(pheno)) {

y <- as.numeric(pheno[,i])
line <- as.factor(pheno$genotype)
block <- as.factor(pheno$Rep)
trt <- as.factor(pheno$treatment)



library(lme4)
yMM <- lmer(y~ (1|line) + (1|block) + (1|trt) + (1|line:block) + (1|line:trt) + (1|block:trt))
summary(yMM)

# estimate BLUPS
yBLUP <- ranef(yMM)
str(yBLUP)

modIntercept <- coef(summary(yMM))[,"Estimate"]
lineEffect <- yBLUP$line

lineBLUP <- modIntercept + lineEffect
#lineBLUP$genotype <- row.names(lineBLUP)
#as.data.frame(lineBLUP)
#colnames(lineBLUP) <- c("BLUP","genotype")

phenoBLUP[[i]] <- lineBLUP

}

phenoBLUPdata <- do.call(cbind,phenoBLUP[6:80])

colnames(phenoBLUPdata) <- c(paste0(colnames(pheno[,6:80]),"_BLUP"))

write.csv(phenoBLUPdata,"C:/Users/natha/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel/traitBLUPs_full.csv")



### RUN LOOP WITH SUBSET BY TREATMENT (treatment==c("control","stress1","stress2"))
## for stress2 treatment, must impute missing values first
phenoBLUP <- list()
phenoControl <- subset(pheno,treatment=="control")
phenoControl$plantno. <- as.factor(phenoControl$plantno.)
phenoControl[6:80] <- sapply(phenoControl[6:80],as.numeric)

library(dplyr)
library(tidyr)
phenoControl_imp <- phenoControl %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) # naive imputation taking mean of column

for (i in 6:ncol(phenoControl_imp)) {
  
  y <- as.numeric(phenoControl_imp[,i])
  line <- as.factor(phenoControl_imp$genotype)
  block <- as.factor(phenoControl_imp$Rep)
  #trt <- as.factor(phenoControl$treatment) # no longer have treatment as an effect due to separation of trials
  
  
  
  library(lme4)
  yMM <- lmer(y~ (1|line) + (1|block))
  summary(yMM)
  
  # estimate BLUPS
  yBLUP <- ranef(yMM)
  str(yBLUP)
  
  modIntercept <- coef(summary(yMM))[,"Estimate"]
  lineEffect <- yBLUP$line
  
  lineBLUP <- modIntercept + lineEffect
  #lineBLUP$genotype <- row.names(lineBLUP)
  #as.data.frame(lineBLUP)
  #colnames(lineBLUP) <- c("BLUP","genotype")
  
  phenoBLUP[[i]] <- lineBLUP
  
}

phenoBLUPdata <- do.call(cbind,phenoBLUP[6:80])

colnames(phenoBLUPdata) <- c(paste0(colnames(pheno[,6:80]),"_BLUP"))

write.csv(phenoBLUPdata,"C:/Users/natha/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel/traitBLUPs_control.csv")



### Make a correlation plot (All Trait BLUPS y-axis VS. All Trait Mean Values x-axis)

controlBLUPS <- read.csv("C:/Users/natha/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel/traitBLUPs_control.csv")
controlBLUPS <- subset(controlBLUPS,environment=="control")


library(dplyr)
library(tidyr)
phenoMEANS <- pheno %>% 
  group_by(treatment,genotype) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

controlMEANS <- subset(phenoMEANS,treatment=="control")
nrow(controlMEANS)
nrow(controlBLUPS)

controlBLUPS <- controlBLUPS[-c(1,3)]
controlMEANS <- controlMEANS[-c(1,3:4)]

BLUPvsMEAN <- merge(controlBLUPS,controlMEANS,by="genotype")
str(BLUPvsMEAN)
BLUPvsMEAN[2:151] <- sapply(BLUPvsMEAN[2:151],as.numeric)
BLUPvsMEANimp<- BLUPvsMEAN %>% 
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))


correlations <- cor(BLUPvsMEANimp[,c(2:151)],method="spearman")
correlations <- correlations[-c(1:75),]
correlations <- correlations[,-c(76:150)]

library(corrplot)
corrplot(correlations,method="color",tl.cex=.45)



#levels <- colnames(BLUPvsMEAN[,c(2:151)])

#library(dplyr)
#library(tidyr)
#library(ggplot2)

#corr_long <- correlations %>% 
#  data.frame() %>% 
#  mutate(row=factor(rownames(.),levels=levels),
#         rowid=as.numeric(row)) %>% 
#  pivot_longer(-c(row,rowid),names_to="col") %>% 
#  mutate(col=factor(col,levels=levels),
#         colid=as.numeric(col))

#ggplot(corr_long, aes(col, row)) +
#  geom_point(aes( fill = value), 
#             data = ~filter(.x, rowid > colid), shape = 21) +
#  geom_text(aes(label = scales::number(value, accuracy = .01), color = abs(value)), 
#            data = ~filter(.x, rowid < colid), size = 8 / .pt) +
#  scale_x_discrete(labels = ~ attr(.x, "pos"), drop = FALSE) +
#  scale_y_discrete(labels = ~ paste0(.x, " (", attr(.x, "pos"), ")"), drop = FALSE) +
#  scale_fill_viridis_c(limits = c(-1, 1)) +
#  scale_color_gradient(low = grey(.8), high = grey(.2)) +
#  coord_equal() +
#  guides(size = "none", color = "none") +
#  theme(legend.position = "bottom", 
#        panel.grid = element_blank(), 
#        axis.ticks = element_blank()) +
#  labs(x = NULL, y = NULL, fill = NULL)

#####################################################################################################
pheno$yield <- as.numeric(pheno$yield)
library(dplyr)
predVSmean <- pheno %>% 
  group_by(genotype) %>% 
  summarize(mean=mean(yield,na.rm=TRUE))
as.data.frame(predVSmean)

predVSmean <- lineBLUP[match(predVSmean$genotype,lineBLUP$genotype)]

predVSmean <- merge(lineBLUP,predVSmean[,c("genotype","mean")],by="genotype")

cor(predVSmean$BLUP,predVSmean$mean) # what is the accuracy of the model?

plot(predVSmean$mean,predVSmean$BLUP,xlim=c(0,15),ylim=c(0,15))
