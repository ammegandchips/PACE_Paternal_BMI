###################################################
# PACE PATERNAL BMI PHASE ONE ANALYSIS PLAN CODE  # 
#                    GEMMA SHARP                  #
#                     29/12/2017                  #
###################################################

####################################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in the PACE Paternal BMI analysis plan.
# If you have insufficient data to complete one or more of the EWAS, you can just skip those models.
# The code also produces .csv files summarising the variables included in the EWASs.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# There are just two inputs required for this analysis:
# 1) pheno: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample(individual) and each column is a different variable. 
#    Necessary variable names are: "pat.bmi", "mat.bmi", pat.active.smoking", "mat.active.smoking", "sex","pat.ses", "pat.age", "mat.age", "parity"
#    If these columns are named differently in your dataset, please rename the columns accordingly
#    Details on how to code these variables are provided in the analysis plan.
# 2) meth: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the sample.id column in pheno.
####################################################################################################################################################################

# Load required packages (if these are not already installed, you will have to install them as the first step)
library(sva)
library(tableone)
library(matrixStats)
library(limma)

# Setup the necessary functions
## Function to remove outliers using the IQR*3 (Tukey) method
IQR.removal <- function(meth.matrix){
  rowIQR <- rowIQRs(meth.matrix, na.rm = T)
  row2575 <- rowQuantiles(meth.matrix, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth.matrix < row2575[,1] - 3 * rowIQR 
  maskU <- meth.matrix > row2575[,2] + 3 * rowIQR 
  meth.matrix[maskL] <- NA
  meth.matrix[maskU] <- NA
  meth.matrix
}
## Function to generate surrogate variables and merge them with the phenotype data (used to adjust for batch)
SVA.generate <- function(meth.matrix, pheno.data, variable.of.interest, model.covariates,n.sv){
  intersecting.samples <- intersect(pheno.data$sample.id,colnames(meth.matrix))
  pheno.data <- na.omit(pheno.data[which(pheno.data$sample.id %in% intersecting.samples),unique(c("sample.id",variable.of.interest,model.covariates))])
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  k = which(is.na(meth.matrix), arr.ind=TRUE)
  meth.matrix[k] = rowMedians(meth.matrix, na.rm=TRUE)[k[,1]]
  mod = model.matrix(reformulate(paste0("pheno.data$",colnames(pheno.data[-1]))))
  mod0 = mod[,-grep(paste0(variable.of.interest,collapse="|"),colnames(mod))]
  sva.ret = sva(meth.matrix, mod=mod, mod0=mod0, n.sv=n.sv)
  SVs = as.data.frame(sva.ret$sv)
  colnames(SVs) <-paste0("sv",1:ncol(SVs))
  cbind(pheno.data,SVs)
}
## Function to run EWAS
ewas.function <-  function(meth.matrix, pheno.data, variable.of.interest){   
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  model.covariates <- colnames(pheno.data)[-which(colnames(pheno.data) %in% c(variable.of.interest,"sample.id"))]
  des = model.matrix(reformulate(paste0("pheno.data$",c(variable.of.interest,model.covariates))))
  fit = lmFit(meth.matrix, des)
  fit.ebayes = eBayes(fit)
  n = rowSums(!is.na(meth.matrix))
  se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$stdev.unscaled))])
  res = data.frame(n=n,
                   coef=fit.ebayes$coefficient[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$coefficient))],
                   se=se,
                   p=fit.ebayes$p.value[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$p.value))])
  res
}

# Set initial parameters
study <- "ALSPAC" #change to your study identifier
timepoint <- "birth" #change depending on the age of the children with methylation samples. Can be "birth", "early_childhood", "late_childhood", "adolescence" or "adult"
cell.names <- if(timepoint=="birth"){
  c("nk","gran","bcell","cd8t","cd4t","mono","nrbc")
}else{
  c("nk","gran","bcell","cd8t","cd4t","mono")
}
traits.and.covariates <- c("pat.bmi", "mat.bmi","sex","pat.ses","pat.age","mat.age","pat.active.smoking","mat.active.smoking","parity")
covariates <- c("pat.age", "pat.active.smoking", "pat.ses", "mat.age", "mat.active.smoking" , "parity", cell.names)

# Load and check phenotype data
pheno <- read.csv("EWAS/pat_bmi/phenofile.alspac.csv",header=TRUE,stringsAsFactors=FALSE) #change filename/location to point to your phenotype file

for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
  print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
               paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
               paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

table(pheno$bio.dad) #checking number of partners that are not biological fathers
pheno <- pheno[which(pheno$bio.dad==1),] #removing partners that are not biological fathers

# Load methylation data (meth)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj") #change to point to your methylation data
# Please perform any cohort-level QC at this stage (e.g. you might want to remove probes with a high detection P-value)

# IQR*3 method to remove outliers (if this has not already been applied to your data)
meth <- IQR.removal(meth)

# Generate surrogate variables for technical batch and merge with pheno data to create the phenotype dataframes for the mutually-adjusted model
pheno.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.bmi","mat.bmi"), model.covariates = c(covariates,"sex"),n.sv=20)

#Create variables for BMI categories, which will be used to summarise data
pheno.mutual$pat.bmi.cat <- NA
pheno.mutual$pat.bmi.cat[pheno.mutual$pat.bmi<=18.5] <-"underweight"
pheno.mutual$pat.bmi.cat[pheno.mutual$pat.bmi>18.5 & pheno.mutual$pat.bmi<=25] <-"normal"
pheno.mutual$pat.bmi.cat[pheno.mutual$pat.bmi>25 & pheno.mutual$pat.bmi<=30] <-"overweight"
pheno.mutual$pat.bmi.cat[pheno.mutual$pat.bmi>30] <-"obese"
pheno.mutual$mat.bmi.cat <- NA
pheno.mutual$mat.bmi.cat[pheno.mutual$mat.bmi<=18.5] <-"underweight"
pheno.mutual$mat.bmi.cat[pheno.mutual$mat.bmi>18.5 & pheno.mutual$mat.bmi<=25] <-"normal"
pheno.mutual$mat.bmi.cat[pheno.mutual$mat.bmi>25 & pheno.mutual$mat.bmi<=30] <-"overweight"
pheno.mutual$mat.bmi.cat[pheno.mutual$mat.bmi>30] <-"obese"

# Create BMI Z-scores
pheno.mutual$Zmat.bmi <- as.numeric(scale(pheno.mutual$mat.bmi,center=TRUE,scale=TRUE))
pheno.mutual$Zpat.bmi <- as.numeric(scale(pheno.mutual$pat.bmi,center=TRUE,scale=TRUE))

#Create the phenotype dataframes for the minimally-adjusted and sex stratified EWASs
pheno.minimal.pat <-pheno.mutual[,-which(colnames(pheno.mutual) %in% c("mat.bmi","Zmat.bmi"))]
pheno.minimal.mat <-pheno.mutual[,-which(colnames(pheno.mutual) %in% c("pat.bmi","Zpat.bmi"))]
pheno.mutual.boys.only <- pheno.mutual[which(pheno.mutual$sex == 0),]
pheno.mutual.girls.only <- pheno.mutual[which(pheno.mutual$sex == 1),]

# Summarise pheno data and save summaries as .csv files
setwd("EWAS/pat_bmi/")
mutual.tableone <- as.data.frame(print(CreateTableOne(data=pheno.mutual[,-1],factorVars=c("pat.bmi.cat","mat.bmi.cat","pat.active.smoking","mat.active.smoking","pat.ses","parity","sex"))),stringsAsFactors=FALSE)
mutual.tableone <- rbind(mutual.tableone,cor.test(pheno.mutual$pat.bmi,pheno.mutual$mat.bmi,method="spearman")$estimate,cor.test(pheno.mutual$pat.bmi,pheno.mutual$mat.bmi,method="spearman")$p.value)
row.names(mutual.tableone)[(nrow(mutual.tableone)-1):nrow(mutual.tableone)] <- c("cortest.rho","cortest.p")
write.csv(mutual.tableone,file=paste0(study,".patbmi.mutual.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.minimal.pat <- ewas.function(meth, pheno.minimal.pat[,!colnames(pheno.minimal.pat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zpat.bmi")
ewas.res.minimal.mat <- ewas.function(meth, pheno.minimal.mat[,!colnames(pheno.minimal.mat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zmat.bmi")
ewas.res.mutual <- ewas.function(meth, pheno.mutual[,!colnames(pheno.mutual) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
ewas.res.mutual.boys.only <- ewas.function(meth, pheno.mutual.boys.only[,!colnames(pheno.mutual.boys.only) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
ewas.res.mutual.girls.only <- ewas.function(meth, pheno.mutual.girls.only[,!colnames(pheno.mutual.girls.only) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
 
# Save EWAS results as an Rdata file
save(list=intersect(ls(),
                    c("ewas.res.minimal.pat",
                      "ewas.res.minimal.mat",
                      "ewas.res.mutual",
                      "ewas.res.mutual.boys.only",
                      "ewas.res.mutual.girls.only")),
     file=paste0(study,".patbmi.ewasresults.",timepoint,".Rdata"))
