# QC of individual cohort EWAS results

#Set wd

setwd("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/")

#Create a key for file/object/results names

key <- data.frame(result=rep(c("min.pat","min.mat","min.patmat","min.matpat",
	"covs.pat","covs.mat","covs.patmat","covs.matpat",
	"boys.patmat","boys.matpat","girls.patmat","girls.matpat"),5),
	model=rep(c("ewas.res.min.pat","ewas.res.min.mat", "ewas.res.min.mutual","ewas.res.min.mutual",
	"ewas.res.covs.pat" , "ewas.res.covs.mat", "ewas.res.covs.mutual","ewas.res.covs.mutual",
	"ewas.res.covs.mutual.boys.only","ewas.res.covs.mutual.boys.only","ewas.res.covs.mutual.girls.only","ewas.res.covs.mutual.girls.only"),5),
	column.suffix=rep(c("pheno.data.Zpat.bmi","pheno.data.Zmat.bmi","pheno.data.Zpat.bmi","pheno.data.Zmat.bmi",
		"pheno.data.Zpat.bmi","pheno.data.Zmat.bmi","pheno.data.Zpat.bmi","pheno.data.Zmat.bmi",
		"pheno.data.Zpat.bmi","pheno.data.Zmat.bmi","pheno.data.Zpat.bmi","pheno.data.Zmat.bmi"),5),
	column.prefix=rep(c("probeid","coef","se","p","n"),each=12),
 	merged.title.coef=c("ewas.res.min.pat","ewas.res.min.mat","ewas.res.min.mutual.coef.pheno.data.Zpat.bmi","ewas.res.min.mutual.coef.pheno.data.Zmat.bmi","ewas.res.covs.pat","ewas.res.covs.mat" ,"ewas.res.covs.mutual.coef.pheno.data.Zpat.bmi","ewas.res.covs.mutual.coef.pheno.data.Zmat.bmi","ewas.res.covs.mutual.boys.only.coef.pheno.data.Zpat.bmi" ,"ewas.res.covs.mutual.boys.only.coef.pheno.data.Zmat.bmi" ,"ewas.res.covs.mutual.girls.only.coef.pheno.data.Zpat.bmi","ewas.res.covs.mutual.girls.only.coef.pheno.data.Zmat.bmi"),
	 merged.title.se=c("ewas.res.min.pat","ewas.res.min.mat","ewas.res.min.mutual.se.pheno.data.Zpat.bmi","ewas.res.min.mutual.se.pheno.data.Zmat.bmi","ewas.res.covs.pat","ewas.res.covs.mat" ,"ewas.res.covs.mutual.se.pheno.data.Zpat.bmi","ewas.res.covs.mutual.se.pheno.data.Zmat.bmi","ewas.res.covs.mutual.boys.only.se.pheno.data.Zpat.bmi" ,"ewas.res.covs.mutual.boys.only.se.pheno.data.Zmat.bmi" ,"ewas.res.covs.mutual.girls.only.se.pheno.data.Zpat.bmi","ewas.res.covs.mutual.girls.only.se.pheno.data.Zmat.bmi"),
	 merged.title.p=c("ewas.res.min.pat","ewas.res.min.mat","ewas.res.min.mutual.p.pheno.data.Zpat.bmi","ewas.res.min.mutual.p.pheno.data.Zmat.bmi","ewas.res.covs.pat","ewas.res.covs.mat" ,"ewas.res.covs.mutual.p.pheno.data.Zpat.bmi","ewas.res.covs.mutual.p.pheno.data.Zmat.bmi","ewas.res.covs.mutual.boys.only.p.pheno.data.Zpat.bmi" ,"ewas.res.covs.mutual.boys.only.p.pheno.data.Zmat.bmi" ,"ewas.res.covs.mutual.girls.only.p.pheno.data.Zpat.bmi","ewas.res.covs.mutual.girls.only.p.pheno.data.Zmat.bmi"),
	stringsAsFactors=FALSE)

#Read in all the EWAS results

extractEWASres <- function(filename){
	load(filename)
	cohort.res <- mget(intersect(unique(key$model),ls()))
	names(cohort.res) <- intersect(unique(key$model),ls())
	cohort.res
}

ALSPAC <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/alspac/results/ALSPAC.patbmi.ewasresults.late_childhood.Rdata")
GenerationR <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/GenR/GenR.patbmi.ewasresults.late_childhood.Rdata")
CHAMACOS <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CHAMACOS_results/CHAMACOS.patbmi.ewasresults.late_childhood.Rdata")
HELIX <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/HELIX/HELIX.patbmi.ewasresults.late_childhood.Rdata")

All.EWAS <- list(ALSPAC,CHAMACOS,GenerationR,HELIX)
names(All.EWAS) <- c("ALSPAC","CHAMACOS","GenerationR","HELIX")

#Add probe ID

addprobeID <- function(ewas.dataframe){
	ewas.dataframe$probeid <- row.names(ewas.dataframe)
	ewas.dataframe
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,addprobeID))

# Remove control probes, remove polymorphic and cross-reactive sites

require(meffil)
annotation <- meffil.get.features("epic") # 867530
annotation <- annotation[-which(annotation$type=="control"),] #remove control probes (485577)
annotation <- annotation[-which(annotation$target=="snp"),] # remove snp probes (485512)
annotation <- annotation[-which(annotation$snp.exclude=="TRUE"),] # remove probes on SNPs (395016)
chen <- read.csv("chen_list.csv",stringsAsFactors=FALSE,header=TRUE)
annotation <- annotation[-which(annotation$name%in%na.omit(chen[,1])),] # remove cross-hybridizing probes (chen) 372993
annotation <- annotation[-which(annotation$chromosome%in%c("chrX","chrY")),] 

filterprobes <- function(ewas.dataframe){
	ewas.dataframe <- ewas.dataframe[intersect(ewas.dataframe$probeid,annotation$name),]
	ewas.dataframe
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,filterprobes))
