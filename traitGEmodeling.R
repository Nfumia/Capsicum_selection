##### Modeling GxE in Genomic Prediction

WorkDir <- "C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel"
setwd(WorkDir)

##Source in functions to be used
source("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Workshops/SummerInstitutePB/gsBootcamp/scripts/bootcamp_functions.R")
gc()



##Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "pepper_core_431.vcf"
genotypes_VCF <- read.table(infileVCF)
vcf <- read.vcfR(infileVCF, verbose = FALSE)
vcf

gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
rownames(gt2) <- rownames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)

gt2d_num<- apply(gt2d,2,as.numeric)
rownames(gt2d_num)<- rownames(gt2d)
geno_num <- t(gt2d_num)
dim(geno_num)
rm(list=grep("gt2",ls(),value=TRUE))


##Filter markers on % missing
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.2)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num

##Filter individuals on % missing
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if(length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2


##Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num2[, -ndx3] else geno_num4 <- geno_num3

dim(geno_num4)


pheno <- read.csv("univariatePHENOBLUP/traitBLUPs_combined.csv")

geno_num4_x <- cbind(rownames(geno_num4),geno_num4)
colnames(geno_num4_x)[1]<- "g2pname"

### Check strain names have same format in pheno and geno 
pheno[,1] <- gsub("[-.]","",pheno[,1])
geno_num4_x[,1] <- gsub("[-.]","",geno_num4_x[,1])

## Merge Geno and Pheno Data
colnames(pheno)[1]<- "g2pname"
Data <- merge(geno_num4_x,pheno,by="g2pname",all=TRUE)

## Remove with missing yield_blup values 

YldNA_Indices <- which(is.na(Data$yield_BLUP))
if(length(YldNA_Indices) >0){Data_Sub <- Data[-YldNA_Indices,]}else{Data_Sub <- Data}


genoStrain <- unique(as.character(geno_num4_x[,"g2pname"]))

genoStrainIndices <- which(Data_Sub[,"g2pname"] %in% genoStrain)
length(genoStrainIndices)
genoIndices <- grep("S",colnames(geno_num4_x))
initGenoIndx <- genoIndices[1]
finalGenoIndx <- genoIndices[length(genoIndices)]
phenoIndices <- c(1,c((finalGenoIndx+1):ncol(Data_Sub)))

pheno_sub <- Data_Sub[genoStrainIndices,phenoIndices]
geno_num4b <- Data_Sub[genoStrainIndices,c(1,genoIndices)]

uniqueStrainIndices<- which(!duplicated(geno_num4b[,"g2pname"]))

if(length(uniqueStrainIndices)>0) {geno_num5 <- geno_num4b[uniqueStrainIndices,]}else{geno_num5 <- geno_num4b}

dim(geno_num5)

rm(geno_num4b)
rm(geno_num4)
rm(geno_num3)
rm(geno_num2)

### set 'yield' colname to 'Yield_blup'

#yldCol <- which(colnames(pheno_sub) %in% "yield")
#colnames(pheno_sub)[yldCol] <- "Yield_blup" 



### Select 3 environs with largest number of evaluations (lines)  

#env_sub <-  names(which(table(pheno_sub[,"environ"])>5100)[1:3])

#env_sub_indices <- which(pheno_sub[,"environ"] %in% env_sub)

## Subset Data and Geno tables 
#DT <- pheno_sub[env_sub_indices,]


DT <- pheno_sub; DT$environment <- as.factor(DT$environment)

dim(DT)

#### Impute genotable using markov function from 'NAM' package 

geno_imp <- markov(apply(geno_num5[,-1],2,as.numeric))
rownames(geno_imp) <- geno_num5[,"g2pname"]
dim(geno_imp)

### 
env_geno_sub_indices <- which(rownames(geno_imp) %in% unique(DT[,"g2pname"]))
geno_imp_sub <- geno_imp[env_geno_sub_indices,]

dim(geno_imp_sub)

K_rr <- A.mat(geno_imp_sub)
colnames(K_rr) <-rownames(geno_imp_sub)
rownames(K_rr) <- rownames(geno_imp_sub)
A <- K_rr
dim(A)


#A_Sub <- A[1:500,1:500]
#DT_Sub <- DT[which(DT[,"strain"] %in% rownames(A_Sub)),]

### Generate Identity Matrix for environments
E <- diag(length(unique(DT$environment)))
rownames(E) <- colnames(E) <- unique(DT$environ)
dim(E)

### Same set of strains in each of the environments 

rmStrains <- names(which(table(DT[,"g2pname"]) <3))
#rmStrains <- names(which(table(DT[,"g2pname"]) >3))
DT_Sub <- DT[-which(DT[,"g2pname"] %in% rmStrains),]

A_Sub <- A[-which(rownames(A) %in% rmStrains),-which(rownames(A) %in% rmStrains)]
dim(A_Sub)

#### Exercise 1 : 


#### 1a) Main Effects Model 

fitMain <- mmer(yield_BLUP~environment-1,
                random=~vsr(g2pname,Gu=A_Sub),
                rcov=~units,
                data=DT_Sub,verbose=FALSE)
summary(fitMain)



m <- model.matrix(~ environment-1,data=DT_Sub)
m_beta <- m %*% as.numeric(fitMain$Beta[,3]) 
PredMain <- m_beta+fitMain$U$`u:g2pname`$yield_BLUP
cor(PredMain,DT_Sub[,"yield_BLUP"]) 
plot(PredMain,DT_Sub[,"yield_BLUP"])

#### 1b) Compound Symmetry Var-Covar 
E <- diag(length(unique(DT_Sub$environment)))
rownames(E) <- colnames(E) <- unique(DT_Sub$environment)

EA <- kronecker(E,A_Sub, make.dimnames = TRUE)
DT_Sub$environment <- as.factor(DT_Sub$environment)
DT_Sub$g2pname <- as.factor(DT_Sub$g2pname)

fitCS <- mmer(biomassfinalplant_BLUP~environment-1,
              random= ~ vsr(g2pname,Gu=A_Sub) + vsr(environment:g2pname,Gu=EA),
              rcov= ~ units,
              data=DT_Sub, verbose = FALSE)
summary(fitCS)


m <- model.matrix(~ environment-1 ,data=DT_Sub)
m_beta <- m %*% as.numeric(fitCS$Beta[,3]) 
PredCS <- m_beta+fitCS$U$`u:environment:g2pname`$biomassfinalplant_BLUP
cor(PredCS,DT_Sub[,"biomassfinalplant_BLUP"]) 

df <- cbind(DT_Sub[,c("environment","biomassfinalplant_BLUP")],PredCS)


colors <- c("black","red","blue")
plot(df$PredCS,df$biomassfinalplant_BLUP,col=colors[df$environment],xlab="GxE Predicted Plant Biomass",ylab="Plant Biomass BLUP",main="Correlation of BLUP and GxE Predictions",sub="Compound Symmetry Variance-Covariance")
abline(lm(df$PredCS~df$biomassfinalplant_BLUP))
legend(-12000,155000,legend=c("Control","High Temp","Severe Temp"),
       col=c("black","red","blue"),pch=1)
text(-1000,110000,expression(r^2==0.746))

















#### 1c) Compound Symmetry Var-Covar with heterogeneous variance


fitCSDG <- mmer(yield_BLUP~environment-1,
                random=~vsr(g2pname,Gu=A_Sub) +vsr(dsr(environment),g2pname,Gu=A_Sub),
                rcov=~units,
                data=DT_Sub,verbose=FALSE) 

summary(fitCSDG)

m2 <- cbind(c(rep(1,nrow(DT_Sub)/3),rep(0,2*nrow(DT_Sub)/3)),c(rep(0,nrow(DT_Sub)/3),rep(1,nrow(DT_Sub)/3),rep(0,nrow(DT_Sub)/3)),
            c(rep(0,nrow(DT_Sub)/3),rep(0,nrow(DT_Sub)/3),rep(1,nrow(DT_Sub)/3)))

m_beta <- m2 %*% as.numeric(fitCSDG$Beta[,3]) 
length(m_beta)
m_env_strain <- do.call(cbind,lapply(fitCSDG$U,function(x) x$yield_BLUP))
dim(m_env_strain)
envStrain_blup <-c(m_env_strain[,2:4])                              

strain_blup <- rep(fitCSDG$U$`u:g2pname`$yield_BLUP,3)
length(strain_blup)

PredCSDG <- m_beta+strain_blup+envStrain_blup

indES <-  sort.int(as.numeric(DT_Sub[,"g2pname"]),decreasing=FALSE,index.return=TRUE)[[2]]

cor(PredCSDG,DT_Sub[indES,"yield_BLUP"]) 
plot(PredCSDG,DT_Sub[indES,"yield_BLUP"])

