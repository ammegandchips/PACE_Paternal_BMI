#Summary of N for each model

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

extractEWASres <- function(filename){
  load(filename)
  cohort.res <- mget(intersect(unique(key$model),ls()))
  names(cohort.res) <- intersect(unique(key$model),ls())
  cohort.res
}

ALSPAC <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/alspac/results/ALSPAC.patbmi.ewasresults.birth.Rdata")
BIB_asian <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/bib/BIB_asian.patbmi.ewasresults.birth.Rdata")
BIB_white <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/bib/BIB_white.patbmi.ewasresults.birth.Rdata")
GenerationR <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/GenR/GenR.patbmi.ewasresults.birth.Rdata")
ProjectViva <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/VIVA_PATBMI/VIVA_PATBMI/Viva.patbmi.ewasresults.birth.Rdata")
CHAMACOS <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CHAMACOS_results/CHAMACOS.patbmi.ewasresults.birth.Rdata")
INMA.nocombat <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/results_PACE_INMA/0years_nocombat/PACE.patbmi.ewasresults.birthnocombat.Rdata")
INMA.combat <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/results_PACE_INMA/0years_combat/PACE.patbmi.ewasresults.birth.Rdata")
RHEA <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/RHEA/RHEA/RHEA.patbmi.ewasresults.birth.Rdata")
ENVIRONAGE <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/ENVIRONAGE/ENVIRONAGE.patbmi.ewasresults.birth.Rdata")
PICCOLIPIU <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/PICCOLIPIU/PICCOLIPIU/PICCOLIPIU.patbmi.ewasresults.birth.Rdata")
GOYA <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/goya/GOYA.patbmi.ewasresults.birth.Rdata")

All.EWAS <- list(ALSPAC,BIB_asian,BIB_white,CHAMACOS,ENVIRONAGE,GenerationR,GOYA,INMA.nocombat,RHEA,PICCOLIPIU,ProjectViva)
names(All.EWAS) <- c("ALSPAC","BIB_asian","BIB_white","CHAMACOS","ENVIRONAGE","GenerationR","GOYA","INMA","RHEA","PICCOLIPIU","ProjectViva")

cutoff=0.1 #calculated previously in QC 5a

remove.outlying.es <- function(ewas.dataframe){
  cbind(ewas.dataframe[,-grep(colnames(ewas.dataframe),pattern="coef")],
    replace(ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")],abs(ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")])>cutoff,NA))
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,remove.outlying.es))

replace.name<-function(X){
  colnames(X)[grep(colnames(X),pattern="ewas")]<-"coef"
  X
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,replace.name))

extract.Ns <- function(ewas.dataframe){
  summary(ewas.dataframe$n)
}

All.Ns<-lapply(All.EWAS,function(x) do.call(rbind,lapply(x,extract.Ns)))
names(All.Ns) <-c("ALSPAC","BIB_asian", "BIB_white","CHAMACOS", "ENVIRONAGE","GenerationR","INMA","PICCOLIPIU", "ProjectViva", "RHEA")
Ns.covs<-do.call(rbind,lapply(All.Ns,function(x) x[which(row.names(x) =="ewas.res.covs.pat"),]))
Ns.covs.mutual.boys.only<-do.call(rbind,lapply(All.Ns,function(x) x[which(row.names(x) =="ewas.res.covs.mutual.boys.only"),]))
Ns.covs.mutual.girls.only<-do.call(rbind,lapply(All.Ns,function(x) x[which(row.names(x) =="ewas.res.covs.mutual.girls.only"),]))

max.covs<-sum(Ns.covs[,6]) #3079
max.boys<-sum(Ns.covs.mutual.boys.only[,6]) #1557
max.girls<-sum(Ns.covs.mutual.girls.only[,6]) #1495
