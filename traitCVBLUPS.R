### Load libraries, functions, and data ----

# Set directory where you have the data files as working dir.. 
WorkDir <- "C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/SideProjects/Pepper_WorldVeg/pepperVarietySel/selectionIndices/genomicSel"
setwd(WorkDir)

# Source in functions to be used. This also loads in all the libraries you need. Change the path to point R towards where you saved this script.
source("C:/Users/natha/OneDrive/Desktop/UH_Manoa/PhD/Workshops/SummerInstitutePB/gsBootcamp/scripts/bootcamp_functions.R")

# Set type of imputation method for use below, either naive or markov
impMethod <- "markov"

# Import phenotypic data (each BLUP dataset from MM - control, stress1, stress2)
pheno <- read.csv("traitBLUPs_control.csv")
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


### Fit some genomic prediction models to the data ---- univariate GBLUP LOOCV through all traits

ndxShuf <- sample(1:dim(geno_imp)[1], dim(geno_imp)[1]) #shuffling genotypes first because RILs are grouped by Parental lineage

pheno_shuf <- pheno2[ndxShuf, ]
geno_imp_shuf <- geno_imp[ndxShuf, ]

cnt <- 1:floor(length(ndxShuf)/nrow(pheno2))

pred_stor <- vector(length=length(ndxShuf))


G <- A.mat(geno_imp)
Stress2GBLUP <- list() #change 'Control' to 'Stress1' and 'Stress2'
Stress2GCorr <- list()
for (i in 4:ncol(pheno2)) {
  
  y <- as.numeric(pheno2[,i])
  G <- G
  df <- pheno2[,c(1,i)]
  
  # Fit an GBLUP model using the rrBLUP package
  gblupModel <- kin.blup(data=df,geno='g2pname',pheno=colnames(df)[2],K=G)
  gblupGebv <- gblupModel$g
  
  
  
  Stress2GBLUP[[i]] <- gblupGebv
  Stress2GCorr[[i]] <- cor(gblupGebv,pheno2[,i])
  
}

#combine the loop lists in to dataframes
Stress2GBLUP.data <- do.call(cbind,Stress2GBLUP[4:78])
Stress2GCorr.data <- do.call(cbind,Stress2GCorr[4:78])

#apply phenotype names to dataframes
colnames(Stress2GBLUP.data) <- c(paste0(colnames(pheno2[,4:78]),"_GBLUP"))
colnames(Stress2GCorr.data) <- c(paste0(colnames(pheno2[,4:78]),"_GCorr"))
Stress2GCorr.data <- as.data.frame(t(Stress2GCorr.data))
Stress2UniG.accuracy <- mean(Stress2GCorr.data[,1],na.rm=TRUE)


#write.csv(Stress2GCorr.data,"traitGCorrelations_stress2.csv")
#write.csv(Stress2GBLUP.data,"traitGBLUPS_stress2.csv")