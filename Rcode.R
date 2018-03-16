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
#    Necessary variable names are: "pat.bmi", "mat.bmi", pat.active.smoking", "mat.active.smoking", "sex","ses", "pat.age", "mat.age", "parity"
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
traits.and.covariates <- c("pat.bmi", "mat.bmi","sex","ses","pat.age","mat.age","pat.active.smoking","mat.active.smoking","parity")
covariates <- c("pat.age", "pat.active.smoking", "ses", "mat.age", "mat.active.smoking" , "parity", cell.names)

# Load and check phenotype data
pheno <- read.csv("EWAS/pat_bmi/phenofile.alspac.csv",header=TRUE,stringsAsFactors=FALSE) #change filename/location to point to your phenotype file

for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
  print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
               paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
               paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

# Load methylation data (meth)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj") #change to point to your methylation data
# Please perform any cohort-level QC at this stage (e.g. you might want to remove probes with a high detection P-value)

# IQR*3 method to remove outliers (if this has not already been applied to your data)
log.iqr <- data.frame(cpgs = row.names(meth),NAs.before.IQR3 = rowSums(is.na(meth)))
meth <- IQR.removal(meth)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(meth))
save(log.iqr, file=paste0(study,".patbmi.logIQR.",timepoint,".Rdata"))

# Generate surrogate variables for technical batch and merge with pheno data to create the phenotype dataframes for the mutually adjusted models
pheno.covs.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.bmi","mat.bmi"), model.covariates = c(covariates,"sex"),n.sv=20)
pheno.min.mutual <- SVA.generate(meth, pheno, variable.of.interest = c("pat.bmi","mat.bmi"), model.covariates = cell.names,n.sv=20)

#Create variables for BMI categories, which will be used to summarise data
pheno.min.mutual$pat.bmi.cat <- NA
pheno.min.mutual$pat.bmi.cat[pheno.min.mutual$pat.bmi<=18.5] <-"underweight"
pheno.min.mutual$pat.bmi.cat[pheno.min.mutual$pat.bmi>18.5 & pheno.min.mutual$pat.bmi<=25] <-"normal"
pheno.min.mutual$pat.bmi.cat[pheno.min.mutual$pat.bmi>25 & pheno.min.mutual$pat.bmi<=30] <-"overweight"
pheno.min.mutual$pat.bmi.cat[pheno.min.mutual$pat.bmi>30] <-"obese"
pheno.min.mutual$mat.bmi.cat <- NA
pheno.min.mutual$mat.bmi.cat[pheno.min.mutual$mat.bmi<=18.5] <-"underweight"
pheno.min.mutual$mat.bmi.cat[pheno.min.mutual$mat.bmi>18.5 & pheno.min.mutual$mat.bmi<=25] <-"normal"
pheno.min.mutual$mat.bmi.cat[pheno.min.mutual$mat.bmi>25 & pheno.min.mutual$mat.bmi<=30] <-"overweight"
pheno.min.mutual$mat.bmi.cat[pheno.min.mutual$mat.bmi>30] <-"obese"

pheno.covs.mutual$pat.bmi.cat <- NA
pheno.covs.mutual$pat.bmi.cat[pheno.covs.mutual$pat.bmi<=18.5] <-"underweight"
pheno.covs.mutual$pat.bmi.cat[pheno.covs.mutual$pat.bmi>18.5 & pheno.covs.mutual$pat.bmi<=25] <-"normal"
pheno.covs.mutual$pat.bmi.cat[pheno.covs.mutual$pat.bmi>25 & pheno.covs.mutual$pat.bmi<=30] <-"overweight"
pheno.covs.mutual$pat.bmi.cat[pheno.covs.mutual$pat.bmi>30] <-"obese"
pheno.covs.mutual$mat.bmi.cat <- NA
pheno.covs.mutual$mat.bmi.cat[pheno.covs.mutual$mat.bmi<=18.5] <-"underweight"
pheno.covs.mutual$mat.bmi.cat[pheno.covs.mutual$mat.bmi>18.5 & pheno.covs.mutual$mat.bmi<=25] <-"normal"
pheno.covs.mutual$mat.bmi.cat[pheno.covs.mutual$mat.bmi>25 & pheno.covs.mutual$mat.bmi<=30] <-"overweight"
pheno.covs.mutual$mat.bmi.cat[pheno.covs.mutual$mat.bmi>30] <-"obese"

# Create BMI Z-scores
pheno.min.mutual$Zmat.bmi <- as.numeric(scale(pheno.min.mutual$mat.bmi,center=TRUE,scale=TRUE))
pheno.min.mutual$Zpat.bmi <- as.numeric(scale(pheno.min.mutual$pat.bmi,center=TRUE,scale=TRUE))

pheno.covs.mutual$Zmat.bmi <- as.numeric(scale(pheno.covs.mutual$mat.bmi,center=TRUE,scale=TRUE))
pheno.covs.mutual$Zpat.bmi <- as.numeric(scale(pheno.covs.mutual$pat.bmi,center=TRUE,scale=TRUE))

#Create the phenotype dataframes for the non-mutually-adjusted and sex stratified EWASs
pheno.min.pat <-pheno.min.mutual[,-which(colnames(pheno.min.mutual) %in% c("mat.bmi","Zmat.bmi","mat.bmi.cat"))]
pheno.min.mat <-pheno.min.mutual[,-which(colnames(pheno.min.mutual) %in% c("pat.bmi","Zpat.bmi","pat.bmi.cat"))]
pheno.covs.pat <-pheno.covs.mutual[,-which(colnames(pheno.covs.mutual) %in% c("mat.bmi","Zmat.bmi"))]
pheno.covs.mat <-pheno.covs.mutual[,-which(colnames(pheno.covs.mutual) %in% c("pat.bmi","Zpat.bmi"))]
pheno.covs.mutual.boys.only <- pheno.covs.mutual[which(pheno.covs.mutual$sex == 0),]
pheno.covs.mutual.girls.only <- pheno.covs.mutual[which(pheno.covs.mutual$sex == 1),]

# Summarise pheno data and save summaries as .csv files
setwd("EWAS/pat_bmi/")
min.mutual.tableone <- as.data.frame(print(CreateTableOne(data=pheno.min.mutual[,-1],factorVars=c("pat.bmi.cat","mat.bmi.cat"))),stringsAsFactors=FALSE)
min.mutual.tableone <- rbind(min.mutual.tableone,suppressWarnings(cor.test(pheno.min.mutual$pat.bmi,pheno.min.mutual$mat.bmi,method="spearman")$estimate),suppressWarnings(cor.test(pheno.min.mutual$pat.bmi,pheno.min.mutual$mat.bmi,method="spearman")$p.value))
covs.mutual.tableone <- as.data.frame(print(CreateTableOne(data=pheno.covs.mutual[,-1],factorVars=c("pat.bmi.cat","mat.bmi.cat","pat.active.smoking","mat.active.smoking","ses","parity","sex"))),stringsAsFactors=FALSE)
covs.mutual.tableone <- rbind(covs.mutual.tableone,suppressWarnings(cor.test(pheno.covs.mutual$pat.bmi,pheno.covs.mutual$mat.bmi,method="spearman")$estimate),suppressWarnings(cor.test(pheno.covs.mutual$pat.bmi,pheno.covs.mutual$mat.bmi,method="spearman")$p.value))
row.names(min.mutual.tableone)[(nrow(min.mutual.tableone)-1):nrow(min.mutual.tableone)] <- c("cortest.rho","cortest.p")
row.names(covs.mutual.tableone)[(nrow(covs.mutual.tableone)-1):nrow(covs.mutual.tableone)] <- c("cortest.rho","cortest.p")
write.csv(min.mutual.tableone,file=paste0(study,".patbmi.min.mutual.summary.",timepoint,".csv"))
write.csv(covs.mutual.tableone,file=paste0(study,".patbmi.covs.mutual.summary.",timepoint,".csv"))

# Test associations between paternal BMI and cell types
cells <- pheno.min.mutual[,which(colnames(pheno.min.mutual) %in% cell.names)]
cells.res <- t(apply(cells,2,function(x) summary(lm(x ~ pheno.min.mutual$pat.bmi))$coef[2,]))
write.csv(cells.res,file=paste0(study,".patbmi.cells.res.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.min.pat <- ewas.function(meth, pheno.min.pat[,!colnames(pheno.min.pat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zpat.bmi")
ewas.res.min.mat <- ewas.function(meth, pheno.min.mat[,!colnames(pheno.min.mat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zmat.bmi")
ewas.res.min.mutual <- ewas.function(meth, pheno.min.mutual[,!colnames(pheno.min.mutual) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
ewas.res.covs.pat <- ewas.function(meth, pheno.covs.pat[,!colnames(pheno.covs.pat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zpat.bmi")
ewas.res.covs.mat <- ewas.function(meth, pheno.covs.mat[,!colnames(pheno.covs.mat) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = "Zmat.bmi")
ewas.res.covs.mutual <- ewas.function(meth, pheno.covs.mutual[,!colnames(pheno.covs.mutual) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
ewas.res.covs.mutual.boys.only <- ewas.function(meth, pheno.covs.mutual.boys.only[,!colnames(pheno.covs.mutual.boys.only) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
ewas.res.covs.mutual.girls.only <- ewas.function(meth, pheno.covs.mutual.girls.only[,!colnames(pheno.covs.mutual.girls.only) %in% c("sex","pat.bmi.cat","mat.bmi.cat", "pat.bmi", "mat.bmi")], variable.of.interest = c("Zpat.bmi","Zmat.bmi"))
 
# Save EWAS results as an Rdata file
save(list=intersect(ls(),
                    c("ewas.res.min.pat",
                      "ewas.res.min.mat",
                      "ewas.res.min.mutual",
                      "ewas.res.covs.pat",
                      "ewas.res.covs.mat",
                      "ewas.res.covs.mutual",
                      "ewas.res.covs.mutual.boys.only",
                      "ewas.res.covs.mutual.girls.only")),
     file=paste0(study,".patbmi.ewasresults.",timepoint,".Rdata"))
