#Read EWAS meta-analysis results and prepare for further analysis

models <-c("min.pat","min.mat","min.patmat","min.matpat","covs.pat",
"covs.mat","covs.patmat","covs.matpat","boys.patmat","boys.matpat",
"girls.patmat","girls.matpat")

setwd("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/")

require(meffil)
annotation <- meffil.get.features("450k")

require(data.table)

read.meta.results <- function(model,tp) {
as.data.frame(fread(paste0("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/",tp,".outliersremoved.",model,"1.txt")))
}

list.of.results<-lapply(models,read.meta.results,tp="birth")

names(list.of.results)<-models

filterprobes <- function(ewas.dataframe){
	ewas.dataframe <- ewas.dataframe[which(ewas.dataframe$MarkerName %in% annotation$name),]
	ewas.dataframe
}

list.of.results <- lapply(list.of.results,filterprobes)
