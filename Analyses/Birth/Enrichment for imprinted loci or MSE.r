#Enrichment for imprinted and metastable epialleles

imprinted <-read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/partida.imprinted.csv",stringsAsFactors=FALSE)[,1]
epialleles <-read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/epiallele cpgs.csv",stringsAsFactors=FALSE)[,1]

myQQ <- function(p.sva, ci = 0.95,Title) {
  n  <- length(p.sva)
  df <- data.frame(
    observed = -log10(p.sva),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
	ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2) +
    ggtitle(Title) +
	theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(log10Pe) +
    ylab(log10Po)  
} 

#imprinted
png("EnrichmentQQ.imprinted.patmatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.patmat$Pvalue[which(list.of.results.het.removed$covs.patmat$MarkerName %in% imprinted)]),Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for imprinted CpGs")
dev.off()
png("EnrichmentQQ.imprinted.matpatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.matpat$Pvalue[which(list.of.results.het.removed$covs.matpat$MarkerName %in% imprinted)]),Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for imprinted CpGs")
dev.off()

#MSE
png("EnrichmentQQ.epialleles.patmatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.patmat$Pvalue[which(list.of.results.het.removed$covs.patmat$MarkerName %in% epialleles)]),Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for metastable epialleles")
dev.off()
png("EnrichmentQQ.epialleles.matpatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.matpat$Pvalue[which(list.of.results.het.removed$covs.matpat$MarkerName %in% epialleles)]),Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for metastable epialleles")
dev.off()

#Summary of P<0.05 in these regions
summary(list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$Pvalue<0.05),"MarkerName"] %in% imprinted)#1
summary(list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$Pvalue<0.05),"MarkerName"] %in% imprinted)#14
summary(list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$Pvalue<0.05),"MarkerName"] %in% epialleles)#8
summary(list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$Pvalue<0.05),"MarkerName"] %in% epialleles)#4

#Denominators
length(intersect(imprinted,list.of.results.het.removed$covs.patmat$MarkerName))#104
length(intersect(epialleles,list.of.results.het.removed$covs.patmat$MarkerName))#110

#Enrichment for imprinted CpGs (now stratified by maternally and paternally derived)

myQQ2 <- function(dat, ci = 0.95,Title) {
require(meffil)
  n  <- nrow(dat)
  dat<-merge(dat,partida,by.x="MarkerName",by.y="CpG",all=FALSE)
  dat$CpG.gene <- paste0(dat$MarkerName," [",dat$Closest.known.imprinted.gene,"]")
  df <- data.frame(
    observed = -log10(sort(dat$Pvalue)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)),
	CpG = dat$CpG.gene
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
	geom_text(aes(expected, observed, label=CpG), angle=90,hjust=0, nudge_y=0.05, fontface="italic",size=3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5,colour="blue") +
    geom_line(aes(expected, cupper), linetype = 2,colour="blue") +
	#ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2,colour="blue") +
    ggtitle(Title) +
	theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(log10Pe) +
    ylab(log10Po)  
} 

partida <-read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/partida.sup.T3.csv",stringsAsFactors=FALSE)
require(data.table)
pat <- as.data.frame(fread("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/geneimprint.paternally.derived.csv",stringsAsFactors=F))
mat <- as.data.frame(fread("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/geneimprint.maternally.derived.csv",stringsAsFactors=F))
pat.regions<-c(pat[,1],unlist(strsplit(pat$Aliases,split=", ")))
mat.regions<-c(mat[,1],unlist(strsplit(mat$Aliases,split=", ")))
partida$parent<-NA
partida$parent[which(partida$Closest.known.imprinted.gene %in% pat.regions)]<-"pat"
partida$parent[which(partida$Closest.known.imprinted.gene %in% mat.regions)]<-"mat"

png("EnrichmentQQ.paternallyimprinted.patmatbmi.png",width=700,height=500)
myQQ2(dat=list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.matpatbmi.png",width=700,height=500)
myQQ2(dat=list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

png("EnrichmentQQ.paternallyimprinted.matpatbmi.png",width=700,height=500)
myQQ2(dat = list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.patmatbmi.png",width=700,height=500)
myQQ2(dat = list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()
