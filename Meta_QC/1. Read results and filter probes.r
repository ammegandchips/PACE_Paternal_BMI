#Read meta-analysis results and filter probes so that only 450k probes remain (BiB used EPIC)

models <-c("min.pat","min.mat","min.patmat","min.matpat","covs.pat",
"covs.mat","covs.patmat","covs.matpat","boys.patmat","boys.matpat",
"girls.patmat","girls.matpat")

tp <- "birth" #or whatever

setwd("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/")

require(meffil)
require(data.table)

annotation <- meffil.get.features("450k")

read.meta.results <- function(model,tp) {
as.data.frame(fread(paste0("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/",tp,".",model,"1.txt")))
}

list.of.results<-lapply(models,read.meta.results,tp="birth")

names(list.of.results)<-models

#Filter just to 450k (BIB had EPIC data)

filterprobes <- function(ewas.dataframe){
	ewas.dataframe <- ewas.dataframe[which(ewas.dataframe$MarkerName %in% annotation$name),]
	ewas.dataframe
}

list.of.results <- lapply(list.of.results,filterprobes)
