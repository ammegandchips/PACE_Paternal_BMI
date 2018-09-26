#Summarise results and hits
#Including number of hits in imprinted regions (from Partida et al.) and P-values<0.05 in imprinted regions

summarise.hits <- function(meta.res){
bonf <- length(meta.res$Pvalue[which(p.adjust(meta.res$Pvalue,method="bonferroni")<0.05)])
bonf.es <- range(abs(meta.res$Effect[which(p.adjust(meta.res$Pvalue,method="bonferroni")<0.05)]))
fdr <- length(meta.res$Pvalue[which(p.adjust(meta.res$Pvalue,method="fdr")<0.05)])
fdr.es <- range(abs(meta.res$Effect[which(p.adjust(meta.res$Pvalue,method="fdr")<0.05)]))
p5 <-  length(meta.res$Pvalue[which(meta.res$Pvalue<1e-5)])
p5.es <-  range(abs(meta.res$Effect[which(meta.res$Pvalue<1e-5)]))
nrow.imprinted <- nrow(meta.res[which(meta.res$MarkerName %in% imprinted[,1]),])
meta.res.imprinted <- meta.res[which(meta.res$MarkerName %in% imprinted[,1]),]
bonf.imprinted <- length(meta.res.imprinted$Pvalue[which(p.adjust(meta.res.imprinted$Pvalue,method="bonferroni")<0.05)])
bonf.es.imprinted <- range(abs(meta.res.imprinted$Effect[which(p.adjust(meta.res.imprinted$Pvalue,method="bonferroni")<0.05)]))
fdr.imprinted <- length(meta.res.imprinted$Pvalue[which(p.adjust(meta.res.imprinted$Pvalue,method="fdr")<0.05)])
fdr.es.imprinted <- range(abs(meta.res.imprinted$Effect[which(p.adjust(meta.res.imprinted$Pvalue,method="fdr")<0.05)]))
p5.imprinted <-  length(meta.res.imprinted$Pvalue[which(meta.res.imprinted$Pvalue<0.05)])
p5.es.imprinted <-  range(abs(meta.res.imprinted$Effect[which(meta.res.imprinted$Pvalue<0.05)]))
list(bonf,bonf.es,fdr,fdr.es,p5,p5.es,
nrow.imprinted,bonf.imprinted,bonf.es.imprinted,fdr.imprinted,fdr.es.imprinted,p5.imprinted,p5.es.imprinted)
}

imprinted <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/partida.imprinted.csv", header=TRUE,stringsAsFactors=FALSE)

Summary.of.hits <-lapply(list.of.results,summarise.hits)
Summary.of.hits <-as.data.frame(do.call(rbind, Summary.of.hits))
colnames(Summary.of.hits) <- c("FE.N.Bonf","FE.Range.ES.Bonf","FE.N.FDR","FE.Range.ES.FDR","FE.N.P<1e-5","FE.Range.ES.P5",
"N imprinted","imprinted.N.Bonf","imprinted.Range.ES.Bonf","imprinted.N.FDR","imprinted.Range.ES.FDR","imprinted.N.P0.05","imprinted.Range.ES.0.05")
Summary.of.hits <- as.data.frame(t(apply(Summary.of.hits,1,unlist)))
Summary.of.hits <-do.call(data.frame,lapply(Summary.of.hits, function(x) replace(x, is.infinite(x),NA)))
row.names(Summary.of.hits) <- names(list.of.results)

time_point<-"birth" #or whatever
write.csv(Summary.of.hits,paste0("meta_models.summary_of_hits.",time_point,".csv"),row.names=TRUE)
