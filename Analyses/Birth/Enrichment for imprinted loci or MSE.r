myQQ2 <- function(dat, ci = 0.95,Title) {
require(meffil)
  lambda <- Lambda(dat$Pvalue)
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
  geom_text(aes(expected, observed, label=CpG), angle=90,hjust=0, nudge_y=0.05, fontface="italic",size=5) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5,colour="blue") +
    geom_line(aes(expected, cupper), linetype = 2,colour="blue") +
  #ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2,colour="blue") +
    ggtitle(paste0(Title," Lambda = ",round(lambda,2))) +
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
myQQ2(dat=list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.matpatbmi.png",width=700,height=500)
myQQ2(dat=list.of.results$covs.matpat[which(list.of.results$covs.matpat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

png("EnrichmentQQ.paternallyimprinted.matpatbmi.png",width=700,height=500)
myQQ2(dat = list.of.results$covs.matpat[which(list.of.results$covs.matpat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.patmatbmi.png",width=700,height=500)
myQQ2(dat = list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

png("EnrichmentQQ.paternallyimprinted.patmatbmi.girls.png",width=700,height=500)
myQQ2(dat=list.of.results$girls.patmat[which(list.of.results$girls.patmat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Female offspring: Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.matpatbmi.girls.png",width=700,height=500)
myQQ2(dat=list.of.results$girls.matpat[which(list.of.results$girls.matpat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Female offspring: Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

png("EnrichmentQQ.paternallyimprinted.patmatbmi.boys.png",width=700,height=500)
myQQ2(dat=list.of.results$boys.patmat[which(list.of.results$boys.patmat$MarkerName %in% partida$CpG[which(partida$parent=="pat")]),],Title="Male offspring: Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.maternallyimprinted.matpatbmi.boys.png",width=700,height=500)
myQQ2(dat=list.of.results$boys.matpat[which(list.of.results$boys.matpat$MarkerName %in% partida$CpG[which(partida$parent=="mat")]),],Title="Male offspring: Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

#everything from partida

png("EnrichmentQQ.allpartida.patmatbmi.png",width=700,height=500)
myQQ(p.sva=sort(list.of.results$covs.patmat$Pvalue[which(list.of.results$covs.patmat$MarkerName %in% partida$CpG)]),Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for paternally-imprinted CpGs")
dev.off()
png("EnrichmentQQ.allpartida.matpatbmi.png",width=700,height=500)
myQQ(p.sva=sort(list.of.results$covs.matpat$Pvalue[which(list.of.results$covs.matpat$MarkerName %in% partida$CpG)]),Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for maternally-imprinted CpGs")
dev.off()

# Previous studies have found associations with lower methylation at MEST, PEG3, NNAT, IGF2

myQQ3 <- function(dat, ci = 0.95,Title) {
require(meffil)
  lambda <- Lambda(dat$Pvalue)
  n  <- nrow(dat)
  df <- data.frame(
    sign.effect = sign(dat$Effect),
    observed = dat$Pvalue)
    df<-df[order(df$observed),]
    df$observed <- -log10(df$observed)
    df$observed <- df$sign.effect * df$observed
    df$expected <- -log10(ppoints(n))
    df$clower  <- -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))
    df$cupper   <- -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5,colour="blue") +
    geom_line(aes(expected, cupper), linetype = 2,colour="blue") +
  #ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2,colour="blue") +
    geom_abline(intercept = 0, slope = -1, alpha = 0.5,colour="blue") +
    geom_line(aes(expected, -cupper), linetype = 2,colour="blue") +
  xlim(0,3) +
    geom_line(aes(expected, -clower), linetype = 2,colour="blue") +
    ggtitle(paste0(Title," Lambda = ",round(lambda,2))) +
  theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(log10Pe) +
    ylab(log10Po)  
} 

IGF2 <- annotation[which(annotation$chromosome=="chr11" & (annotation$position>2150342 & annotation$position<2170833)),]
MEST <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>130126012 & annotation$position<130146133)),]
PEG3 <- annotation[which(annotation$chromosome=="chr19" & (annotation$position>57321445 & annotation$position<57352096)),]
NNAT <- annotation[which(annotation$chromosome=="chr20" & (annotation$position>36149607 & annotation$position<36152092)),]

myQQ3(dat=list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),],Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for CpGs at IGF2")

igf2<- df[which(df$MarkerName %in% IGF2$name),]
mest<- df[which(df$MarkerName %in% MEST$name),]
peg3<- df[which(df$MarkerName %in% PEG3$name),]
nnat<- df[which(df$MarkerName %in% NNAT$name),]
dat <- data.frame(region=c(rep("IGF2",nrow(igf2)),
  rep("MEST",nrow(mest)),
  rep("PEG3",nrow(peg3)),
  rep("NNAT",nrow(nnat))),
  direction=c(sign(igf2$Effect),sign(mest$Effect),sign(peg3$Effect),sign(nnat$Effect)),
  significant=c(igf2$Pvalue<0.05,mest$Pvalue<0.05,peg3$Pvalue<0.05,nnat$Pvalue<0.05))
dat$direction.significant <- NA
dat$direction.significant[which(dat$direction=="1" & dat$significant=="FALSE")]<-"hypermethylated P>0.05"
dat$direction.significant[which(dat$direction=="1" & dat$significant=="TRUE")]<-"hypermethylated P<0.05"
dat$direction.significant[which(dat$direction=="-1" & dat$significant=="FALSE")]<-"hypomethylated P>0.05"
dat$direction.significant[which(dat$direction=="-1" & dat$significant=="TRUE")]<-"hypomethylated P<0.05"
P<-ggplot(dat,aes(region))+
geom_bar(aes(fill=as.factor(direction.significant)))+
theme_classic()+
scale_fill_manual(values=c("firebrick3","firebrick4","dodgerblue2","dodgerblue4"))+
theme(plot.title = element_text(hjust = 0.5))+
theme(legend.title=element_blank())+
xlab("Imprinted region")+ylab("Number of CpGs")+
ggtitle("Associations between paternal BMI (with adjustment for maternal BMI)\nand methylation at paternally-derived imprinted regions")

png("Soubry.paternallyimprinted.patmatbmi.png",width=600,height=500)
P
dev.off()
