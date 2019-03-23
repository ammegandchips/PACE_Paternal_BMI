#meta-regression at top sites
#top sites from the childhood EWAS

tp <- "childhood"

models <-c("covs.pat",
"covs.mat","covs.patmat","covs.matpat","boys.patmat","boys.matpat",
"girls.patmat","girls.matpat")

setwd("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/")

require(meffil)
annotation <- meffil.get.features("450k")

require(data.table)

read.meta.results <- function(model,tp) {
as.data.frame(fread(paste0("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/meta/meta_results/",tp,".",model,"1.txt")))
}

list.of.results<-lapply(models,read.meta.results,tp="childhood")

names(list.of.results)<-models

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

ALSPAC <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/alspac/results/ALSPAC.patbmi.ewasresults.late_childhood.Rdata")
GenerationR <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/GenR/GenR.patbmi.ewasresults.late_childhood.Rdata")
CHAMACOS <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CHAMACOS_results/CHAMACOS.patbmi.ewasresults.late_childhood.Rdata")
HELIX <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/HELIX/HELIX.patbmi.ewasresults.late_childhood.Rdata")
INMA <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/results_PACE_INMA/4years_nocombat/PACE.patbmi.ewasresults.early_childhoodnocombat.Rdata")
ProjectViva <- extractEWASres("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/VIVA_PATBMI/VIVA_PATBMI/Viva.patbmi.ewasresults.early_childhood.Rdata")

All.EWAS <- list(ALSPAC,CHAMACOS,GenerationR,HELIX,INMA,ProjectViva)
names(All.EWAS) <- c("ALSPAC","CHAMACOS","GenerationR","HELIX","INMA","ProjectViva")

addprobeID <- function(ewas.dataframe){
  ewas.dataframe$probeid <- row.names(ewas.dataframe)
  ewas.dataframe
}

All.EWAS.top <- lapply(All.EWAS,function(x) lapply(x,addprobeID))

top.sites.covs.pat <- list.of.results$covs.pat[which(list.of.results$covs.pat$Pvalue<1e-5),"MarkerName"]
top.sites.covs.patmat <- list.of.results$covs.patmat[which(list.of.results$covs.patmat$Pvalue<1e-5),"MarkerName"]
top.sites.boys <- list.of.results$boys.patmat[which(list.of.results$boys.patmat$Pvalue<1e-5),"MarkerName"]
top.sites.girls <- list.of.results$girls.patmat[which(list.of.results$girls.patmat$Pvalue<1e-5),"MarkerName"]

All.EWAS.top.covs.pat<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.covs.pat),]))
All.EWAS.top.covs.patmat<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.covs.patmat),]))
All.EWAS.top.boys<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.boys),]))
All.EWAS.top.girls<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.girls),]))

covs.pat<-lapply(All.EWAS.top.covs.pat,function(x) x$ewas.res.covs.mutual)
covs.patmat<-lapply(All.EWAS.top.covs.patmat,function(x) x$ewas.res.covs.mutual)
boys<-lapply(All.EWAS.top.boys,function(x) x$ewas.res.covs.mutual.boys.only)
girls<-lapply(All.EWAS.top.girls,function(x) x$ewas.res.covs.mutual.girls.only)

metareg.covariates<-data.frame(study=names(All.EWAS.top.covs.pat),age=c(7,9,6,8,4,7),in.cord=c(1,1,1,0,1,1))

covs.pat<-do.call(rbind,covs.pat)
covs.patmat<-do.call(rbind,covs.patmat)
boys<-do.call(rbind,boys)
girls<-do.call(rbind,girls)

covs.pat$study <- unlist(lapply(strsplit(row.names(covs.pat),split=".",fixed=TRUE),"[",1))
covs.patmat$study <- unlist(lapply(strsplit(row.names(covs.patmat),split=".",fixed=TRUE),"[",1))
boys$study <- unlist(lapply(strsplit(row.names(boys),split=".",fixed=TRUE),"[",1))
girls$study <- unlist(lapply(strsplit(row.names(girls),split=".",fixed=TRUE),"[",1))

covs.pat <- merge(covs.pat,metareg.covariates,by="study",all=TRUE)
covs.patmat <- merge(covs.patmat,metareg.covariates,by="study",all=TRUE)
boys <- merge(boys,metareg.covariates,by="study",all=TRUE)
girls <- merge(girls,metareg.covariates,by="study",all=TRUE)


mixed.effects.meta.analysis <- function(data,coef.name,se.name){
                              require(metafor)
                              data <- data[,c("probeid","study",coef.name,se.name,"age","in.cord")]
                              colnames(data) <- c("probeid","study","coef","se","age","in.cord")
                data$coef <- data$coef*100
                data$se <- data$se*100
                              res = split(data, f=data$probeid)
                              res = lapply(res, function(x) rma.uni(slab=x$study,yi=x$coef,sei=x$se,mods= ~ x$age,weighted=TRUE))
                              res
                              }
                
covs.pat.res<-mixed.effects.meta.analysis(data=covs.pat,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
covs.patmat.res<-mixed.effects.meta.analysis(data=covs.patmat,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
boys.res<-mixed.effects.meta.analysis(data=boys,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
girls.res<-mixed.effects.meta.analysis(data=girls,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")


meta.regression.plot <- function(X,Y,Title,i){
### adjust margins so the space is better used
par(mar=c(5,5,5,2))
 ### calculate predicted risk ratios for 3 to 9 years
preds <- predict.rma(X[[i]], newmods=c(3:9),intercept=FALSE)
### calculate point sizes by rescaling the standard errors
wi    <- 1/(sqrt(Y$se.pheno.data.Zpat.bmi[which(Y$probeid==unique(Y$probeid)[1])]))
size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))
### plot the risk ratios against absolute latitude
plot(Y$age[which(Y$probeid==unique(Y$probeid)[i])], Y$coef.pheno.data.Zpat.bmi[which(Y$probeid==unique(Y$probeid)[i])], pch=19, cex=size, 
     xlab="Average age of children in cohort", ylab="Change in % methylation per 1SD increase in paternal BMI", ylim=c(min(Y$coef.pheno.data.Zpat.bmi[which(Y$probeid==unique(Y$probeid)[i])],preds$ci.lb)-0.1,max(Y$coef.pheno.data.Zpat.bmi[which(Y$probeid==unique(Y$probeid)[i])],preds$ci.ub)+0.1),
     las=1, bty="l",main=paste0(Title,"\n",unique(Y$probeid)[i],"\nlinear relation between effect estimate and age: B=",round(X[[i]]$b[2],3),"(",round(X[[i]]$ci.lb[2],3),", ",round(X[[i]]$ci.ub[2],3),"), P=",formatC(X[[i]]$pval[2],format="e",digits=2))
)
### add predicted values (and corresponding CI bounds)
lines(3:9, preds$pred,col="blue")
lines(3:9, preds$ci.lb, lty="dashed",col="blue")
lines(3:9, preds$ci.ub, lty="dashed",col="blue")
### dotted line at RR=1 (no difference between groups)
abline(h=0, lty="dotted")
### labels some points in the plot
text(Y$age[which(Y$probeid==unique(Y$probeid)[i])], (Y$coef.pheno.data.Zpat.bmi[which(Y$probeid==unique(Y$probeid)[i])])-0.01, Y$study[which(Y$probeid==unique(Y$probeid)[i])], cex=0.9,pos=1)
}

pdf("metaregression.pat.fromchildhood.pdf",width=13,height=8)
lapply(1:length(covs.pat.res), function(i) meta.regression.plot(X=covs.pat.res, Y=covs.pat, Title="Paternal BMI (without adjustment for maternal BMI)",i=i))
dev.off()

pdf("metaregression.patmat.fromchildhood.pdf",width=13,height=8)
lapply(1:length(covs.patmat.res), function(i) meta.regression.plot(X=covs.patmat.res, Y=covs.patmat, Title="Paternal BMI (with adjustment for maternal BMI)",i=i))
dev.off()

pdf("metaregression.girls.fromchildhood.pdf",width=13,height=8)
lapply(1:length(girls.res), function(i) meta.regression.plot(X=girls.res, Y=girls, Title="Paternal BMI (with adjustment for maternal BMI): Females only",i=i))
dev.off()

pdf("metaregression.boys.fromchildhood.pdf",width=13,height=8)
lapply(1:length(boys.res), function(i) meta.regression.plot(X=boys.res, Y=boys, Title="Paternal BMI (with adjustment for maternal BMI): Males only",i=i))
dev.off()

covs.pat.metareg = do.call(rbind,lapply(covs.pat.res, function(x)unlist(c(x$b[2,],lapply(x[c("se","zval","pval","ci.lb","ci.ub")],"[",2)))))
covs.patmat.metareg = do.call(rbind,lapply(covs.patmat.res, function(x)unlist(c(x$b[2,],lapply(x[c("se","zval","pval","ci.lb","ci.ub")],"[",2)))))
girls.metareg = do.call(rbind,lapply(girls.res, function(x)unlist(c(x$b[2,],lapply(x[c("se","zval","pval","ci.lb","ci.ub")],"[",2)))))
boys.metareg = do.call(rbind,lapply(boys.res, function(x)unlist(c(x$b[2,],lapply(x[c("se","zval","pval","ci.lb","ci.ub")],"[",2)))))

write.csv(covs.pat.metareg,"metaregression.pat.fromchildhood.csv",row.names=TRUE)
write.csv(covs.patmat.metareg,"metaregression.patmat.fromchildhood.csv",row.names=TRUE)
write.csv(girls.metareg,"metaregression.girls.fromchildhood.csv",row.names=TRUE)
write.csv(boys.metareg,"metaregression.boys.fromchildhood.csv",row.names=TRUE)
