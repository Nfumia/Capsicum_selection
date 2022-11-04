### Load libraries, functions, and data ----

# Set directory where you have the data files as working dir.. 
WorkDir <- "C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel"
setwd(WorkDir)

# Source in functions to be used. This also loads in all the libraries you need. Change the path to point R towards where you saved this script.
source("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Workshops/SummerInstitutePB/gsBootcamp/scripts/bootcamp_functions.R")

# Set type of imputation method for use below, either naive or markov
impMethod <- "markov"

# Import phenotypic data (each BLUP dataset from MM - control, stress1, stress2)
pheno <- read.csv("univariatePHENOBLUP/traitBLUPs_stress2.csv")


# Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "pepper_core_431.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)


### Format, manipulate and filter genotype data ----

# Converting VCF file format to numerical matrix format that can be fit in statistical models
gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
gt2d_num <- apply(gt2d, 2, as.numeric)
#Adding row names back in
rownames(gt2d_num) <- rownames(gt2d)
geno_num <- t(gt2d_num)


# Filtering out markers on proportion of missing data
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.5)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num


# Filtering out individuals on proportion of missing data
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if (length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2 


# Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num3[, -ndx3] else geno_num4 <- geno_num3


# Match phenotypic data to marker data
#pheno <- pheno[-c(291,292),] # reduce number of genotypes to ensure k-fold cross validation can work (needs to be multiple of 10)

ndx4 <- match(rownames(geno_num4), pheno$g2pname)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

pheno2 <- pheno[ndx5, ] #ensure length()/nrow() is a multiple of 10
g2pname <- rownames(geno_num4)[ndxNA] #extract genotype names of varieties missing phenotypic data
phenoCand <- as.data.frame(g2pname)
phenoCand[,c(2:78)] <- "NA"
phenoCand[,3] <- "stress2"
colnames(phenoCand) <- colnames(pheno)
pheno <- rbind(pheno2,phenoCand)

pheno <- pheno[-431,]
ndx4 <- match(rownames(geno_num4), pheno$g2pname)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

geno_num5 <- geno_num4[-ndxNA,]
pheno[4:78] <- sapply(pheno[4:78],as.numeric)


# create candidate genotype lists
#ndxCand <- match(rownames(geno_num4),rownames(geno_num5)) #find which genotypes are for candidate population
#ndxCandNA <- which(is.na(ndxCand))
#geno_Cand <- geno_num4[ndxCandNA,]



# Impute genotype data using Markov chain implemented in the NAM package
if (impMethod == "markov") geno_imp <- markov(apply(geno_num5[, -1], 2, as.numeric))
if (impMethod == "markov") rownames(geno_imp) <- rownames(geno_num5)
# Repeat for candidate varieties
#if (impMethod == "markov") geno_Candimp <- markov(apply(geno_Cand[, -1], 2, as.numeric))
#if (impMethod == "markov") rownames(geno_Candimp) <- rownames(geno_Cand)



### Now we use cross-validation to predict candidate phenotypes to the univariate models to account for model overfitting
# Here we are conducting k-fold cross-validation (k=10 so # of genotypes needs to be divisible by 10)
#ndxShuf <- sample(1:dim(geno_imp)[1], dim(geno_imp)[1]) #shuffling genotypes first because RILs are grouped by Parental lineage

#pheno_shuf <- pheno2[ndxShuf, ] #match phenotype rows to ndxShuf
#geno_imp_shuf <- geno_imp[ndxShuf, ] #match genotype rows to ndxShuf

#pred_stor <- vector(length=length(ndxShuf)) #create dummy vector for predicted phenotype of length of genotypes

#G <- A.mat(geno_imp) #generate genomic relationship matrix
#Stress2GBLUP.kcv.CandPred <- list() #change 'Control' to 'Stress1' and 'Stress2'
#Stress2GCorr.kcv.CandPred <- list()
#for (i in 4:ncol(pheno)) {
  
    
  #y <- as.numeric(pheno_shuf[,i])
  #G <- G
#  df <- pheno[,c(1,i)]
  
  # Fit an GBLUP model using the rrBLUP package
#  gblupModel <- kin.blup(data=df,geno='g2pname',pheno=colnames(df)[2],K=G)
#  gblupGebv <- gblupModel$pred #extract GEBV
    
  
#  Stress2GBLUP.kcv.CandPred[[i]] <- gblupGebv
#  Stress2GCorr.kcv.CandPred[[i]] <- cor(gblupGebv,pheno[,i]) #store prediction accuracy by univariate model
  
#}


idx <- match(pheno$g2pname,rownames(geno_imp))
geno_imp <- geno_imp[idx,]

paTrnStor <- list()
phenoTrgtPred <- list()
for (i in 4:ncol(pheno)){
  
  phenoTrn <- pheno[c(1:291), ]
  genoTrn <- geno_imp[c(1:291), ]
  genoTrgt <- geno_imp[c(292:430),]
  phenoTrgt <- pheno[c(292:430),]
  
  ##Train RR-BLUP model using randomly selected training set and predict target population.
  rrModel <- mixed.solve(y=phenoTrn[,i], Z=genoTrn)
  
  mrkEffs <- rrModel$u
  
  ##Use marker effects to calculate genomic estimated breeding values of individuals in training set. Here we are extracting the intercept and adding it back on.
  int <- as.numeric(rrModel$beta)
  rrGebvTrn <- int + genoTrn%*%mrkEffs
  rrGebvTrgt <- int + genoTrgt%*%mrkEffs
  
  paTrnStor[[i]] <- cor(rrGebvTrn, phenoTrn[c(1:291),i])
  phenoTrgtPred[[i]] <- rrGebvTrgt
  
}

#combine the loop lists in to dataframes
Stress2.RRblup.TrgtPredPheno.data <- do.call(cbind,phenoTrgtPred[4:78])
Stress2.RRblup.TrnPredAcc.data <- do.call(cbind,paTrnStor[4:78])

#apply phenotype names to dataframes
colnames(Stress2.RRblup.TrgtPredPheno.data) <- c(paste0(colnames(pheno[,4:78]),"_RRTrgtPred"))
colnames(Stress2.RRblup.TrnPredAcc.data) <- c(paste0(colnames(pheno[,4:78]),"_RRTrnCorr"))
rownames(Stress2.RRblup.TrgtPredPheno.data) <- c(paste0(rownames(genoTrgt)))
Stress2.RRblup.TrnPredAcc.data <- as.data.frame(t(Stress2.RRblup.TrnPredAcc.data))
Stress2.RRblup.TrnPredAcc.data.accuracy <- mean(Stress2.RRblup.TrnPredAcc.data[,1],na.rm=TRUE)

#write.csv(Stress2GCorr.kcv.CandPred.data,"univariateGBLUP/traitGCorrelations_KCV_CandPred_stress2.csv")
#write.csv(Stress2GBLUP.kcv.CandPred.data,"univariateGBLUP/traitGBLUPS_KCV_CandPred_stress2.csv")
