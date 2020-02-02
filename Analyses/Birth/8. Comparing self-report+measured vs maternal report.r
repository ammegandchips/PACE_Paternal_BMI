####### Read in metal results files for self-report/measured and maternal-report, and filter probes as for the main analysis

models <-c("covs.pat","covs.pat.measured.selfreport","covs.pat.maternalreport")

#for childhood, "min" models are omitted (because we didn't use them in the birth analyses. So "models" is set up as:
#models <-c("covs.pat","covs.mat","covs.patmat","covs.matpat","boys.patmat","boys.matpat","girls.patmat","girls.matpat")

tp <- "birth" #or tp <- "childhood"

setwd("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/")

require(meffil)
require(data.table)

annotation <- meffil.get.features("450k")

read.meta.results <- function(model,tp) {
as.data.frame(fread(paste0("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/",tp,".",model,"1.txt")))
}

list.of.results<-lapply(models,read.meta.results,tp=tp)

names(list.of.results)<-models

#Filter just to 450k (BIB had EPIC data)

filterprobes <- function(ewas.dataframe){
	ewas.dataframe <- ewas.dataframe[which(ewas.dataframe$MarkerName %in% annotation$name),]
	ewas.dataframe
}

list.of.results <- lapply(list.of.results,filterprobes)

####### Calculate the number of CpGs with fdr<0.05 and p<1e-5

lapply(list.of.results, function(x) summary(p.adjust(x$Pvalue,method="fdr")<0.05))
#$covs.pat.measured.selfreport
#0
#$covs.pat.maternalreport
#0

lapply(list.of.results, function(x) summary(x$Pvalue<1e-5))
#$covs.pat.measured.selfreport
#   Mode   FALSE    TRUE 
#logical  363245       8 
#$covs.pat.maternalreport
#   Mode   FALSE    TRUE 
#logical  363247       6 

####### Look at concordance between the CpGs with p<1e-5 (but this will be influenced by power)
intersect(list.of.results$covs.pat.measured.selfreport$MarkerName[list.of.results$covs.pat.measured.selfreport$Pvalue<1e-5],
list.of.results$covs.pat.maternalreport$MarkerName[list.of.results$covs.pat.maternalreport$Pvalue<1e-5])
#0
intersect(list.of.results$covs.pat$MarkerName[list.of.results$covs.pat.measured.selfreport$Pvalue<1e-5],
list.of.results$covs.pat.maternalreport$MarkerName[list.of.results$covs.pat.maternalreport$Pvalue<1e-5])
#0
intersect(list.of.results$covs.pat$MarkerName[list.of.results$covs.pat.measured.selfreport$Pvalue<1e-5],
list.of.results$covs.pat.measured.selfreport$MarkerName[list.of.results$covs.pat.maternalreport$Pvalue<1e-5])
#0

####### Calculate correlation between effect estimates at top CpGs from main analysis

cor(list.of.results$covs.pat.measured.selfreport$Effect[list.of.results$covs.pat$Pvalue<1e-5],
list.of.results$covs.pat.maternalreport$Effect[list.of.results$covs.pat$Pvalue<1e-5])
#0.83
       
cor(list.of.results$covs.pat.measured.selfreport$Effect,
list.of.results$covs.pat$Effect)
#0.79

