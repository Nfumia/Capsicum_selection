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
pheno[4:78] <- sapply(pheno[4:78],as.numeric)

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
ndx4 <- match(rownames(geno_num4), pheno$g2pname)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

pheno2 <- pheno[ndx5, ]
geno_num5 <- geno_num4[-ndxNA, ]

# Impute genotype data using Markov chain implemented in the NAM package
if (impMethod == "markov") geno_imp <- markov(apply(geno_num5[, -1], 2, as.numeric))
if (impMethod == "markov") rownames(geno_imp) <- rownames(geno_num5)


### Fit some bayesian genomic prediction models to the data ---- univariate bayes blup through all traits

G <- A.mat(geno_imp)
Stress2BayesCBLUP <- list() #change 'Control' to 'Stress1' and 'Stress2'
Stress2BayesCCorr <- list()
for (i in 4:ncol(pheno2)) {
  
  y <- as.numeric(pheno2[,i])
  
  # Fit an Bayesian model using the BGLR package
  ETA <- list(list(K=NULL, X=geno_imp, model='BayesC', probIn=.10)) 
  model_bglr <- BGLR(y=y, ETA=ETA, burnIn=500, nIter=2000, verbose=FALSE)
  gebv_bglr <- model_bglr$yHat
  
  Stress2BayesCBLUP[[i]] <- gebv_bglr
  Stress2BayesCCorr[[i]] <- cor(gebv_bglr,pheno2[,i])
  
}

#combine the loop lists in to dataframes
Stress2BayesCBLUP.data <- do.call(cbind,Stress2BayesCBLUP[4:78])
Stress2BayesCCorr.data <- do.call(cbind,Stress2BayesCCorr[4:78])

#apply phenotype names to dataframes
colnames(Stress2BayesCBLUP.data) <- c(paste0(colnames(pheno2[,4:78]),"_BayesCBLUP"))
colnames(Stress2BayesCCorr.data) <- c(paste0(colnames(pheno2[,4:78]),"_BayesCCorr"))
Stress2BayesCCorr.data <- as.data.frame(t(Stress2BayesCCorr.data))
Stress2BayesCUniG.accuracy <- mean(Stress2BayesCCorr.data[,1],na.rm=TRUE)


#write.csv(Stress1BayesCCorr.data,"univariateBayesBLUP/traitBayesCCorrelations_stress1.csv")
#write.csv(Stress1BayesCBLUP.data,"univariateBayesBLUP/traitBayesCBLUPS_stress1.csv")


### Now introduce cross-validation to the univariate models to account for model overfitting

ndxShuf <- sample(1:dim(geno_imp)[1], dim(geno_imp)[1]) #shuffling genotypes first because RILs are grouped by Parental lineage

pheno_shuf <- pheno2[ndxShuf, ] #match phenotype rows to ndxShuf
geno_imp_shuf <- geno_imp[ndxShuf, ] #match genotype rows to ndxShuf

pred_stor <- vector(length=length(ndxShuf)) #create dummy vector for predicted phenotype of length of genotypes

G <- A.mat(geno_imp_shuf) #generate genomic relationship matrix
Stress2BayesCBLUP.loocv <- list() #change 'Control' to 'Stress1' and 'Stress2'
Stress2BayesCCorr.loocv <- list()
for (i in 4:ncol(pheno_shuf)) {
  
  cnt <- 1:floor(length(ndxShuf)/291)  
  #G <- G
  
  for (j in 1:291){
    
    y <- as.numeric(pheno_shuf[,i])
    y[cnt] <- NA #hide a phenotypic value at a time
    
    # Fit an Bayesian model using the BGLR package
    ETA <- list(list(K=NULL, X=geno_imp_shuf, model='BayesC', probIn=.10)) 
    model_bglr <- BGLR(y=y, ETA=ETA, burnIn=500, nIter=2000, verbose=FALSE)
    gebv_bglr <- model_bglr$yHat #extract GEBV
    
    pred_stor[cnt] <- gebv_bglr[cnt] #store the LOOCV results
    
    cnt <- cnt + floor(length(ndxShuf)/291) #change cnt to next row (genotype's phenotype) to be excluded from analysis
  }
  
  Stress2BayesCBLUP.loocv[[i]] <- pred_stor
  Stress2BayesCCorr.loocv[[i]] <- cor(pred_stor,pheno_shuf[,i]) #store prediction accuracy by univariate model
  
}

#combine the loop lists in to dataframes
Stress2BayesCBLUP.loocv.data <- do.call(cbind,Stress2BayesCBLUP.loocv[4:78])
Stress2BayesCCorr.loocv.data <- do.call(cbind,Stress2BayesCCorr.loocv[4:78])

#apply phenotype names to dataframes
colnames(Stress2BayesCBLUP.loocv.data) <- c(paste0(colnames(pheno_shuf[,4:78]),"_BayesBLUPLOOCV"))
colnames(Stress2BayesCCorr.loocv.data) <- c(paste0(colnames(pheno_shuf[,4:78]),"_BayesCorrLOOCV"))
rownames(Stress2BayesCBLUP.loocv.data) <- c(paste0(rownames(geno_imp_shuf)))
Stress2BayesCCorr.loocv.data <- as.data.frame(t(Stress2BayesCCorr.loocv.data))
Stress2BayesCCorr.loocv.data.accuracy <- mean(Stress2BayesCCorr.loocv.data[,1],na.rm=TRUE)

write.csv(Stress2BayesCCorr.loocv.data,"univariateBayesBLUP/traitBayesCCorrelations_LOOCV_stress2.csv")
write.csv(Stress2BayesCBLUP.loocv.data,"univariateBayesBLUP/traitBayesCBLUPS_LOOCV_stress2.csv")
