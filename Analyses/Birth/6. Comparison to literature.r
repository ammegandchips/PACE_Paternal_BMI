#In the literature, we found evidence for association with paternal BMI at imprinted regions and at >9000 CpGs in a paper by Donkin et al.
#We also compare to results from our 2018 PACE Maternal BMI paper

#Imprinted regions

## heatmaps

require(meffil)

IGF2 <- annotation[which(annotation$chromosome=="chr11" & (annotation$position>2150342 & annotation$position<2170833)),]
MEST <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>130126012 & annotation$position<130146133)),]
PEG3 <- annotation[which(annotation$chromosome=="chr19" & (annotation$position>57321445 & annotation$position<57352096)),]
NNAT <- annotation[which(annotation$chromosome=="chr20" & (annotation$position>36149607 & annotation$position<36152092)),]
NDN <- annotation[which(annotation$chromosome=="chr15" & (annotation$position>23930554 & annotation$position<23932450)),]
SNRPN <- annotation[which(annotation$chromosome=="chr15" & (annotation$position>25068794 & annotation$position<25664609)),]
SGCE <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>94214536 & annotation$position<94285521)),]
PEG10 <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>94285637 & annotation$position<94299007)),]
MEG3 <- annotation[which(annotation$chromosome=="chr14" & (annotation$position>101245747 & annotation$position<101327368)),]
H19 <- annotation[which(annotation$chromosome=="chr11" & (annotation$position>2016406 & annotation$position<2022700)),]
PLAGL1 <- annotation[which(annotation$chromosome=="chr6" & (annotation$position>144261437 & annotation$position<144385735)),]
GRB10 <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>50657760 & annotation$position<50861159)),]

imprinted.regions <-  rbind(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),])

imprinted.regions$imprinted.region <- c(rep("IGF2 (n probes = 97)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),])),
rep("MEST (n probes = 64)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),])),
rep("PEG3 (n probes = 25)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),])),
rep("NNAT (n probes = 9)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),])),
rep("NDN (n probes = 6)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),])),
rep("SNRPN (n probes = 267)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),])),
rep("SGCE (n probes = 40)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),])),
rep("PEG10 (n probes = 70)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),])),
rep("MEG3 (n probes = 51)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),])),
rep("H19 (n probes = 52)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),])),
rep("PLAGL1 (n probes = 36)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),])),
rep("GRB10 (n probes = 52)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),]))
)

imprinted.regions <- merge(imprinted.regions,annotation,by.x="MarkerName",by.y="name",all.y=F)
imprinted.regions <- imprinted.regions[order(imprinted.regions$position),]


imprinted.regions$ID<-NA
imprinted.regions$ID[imprinted.regions$imprinted.region=="IGF2 (n probes = 97)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% IGF2$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="MEST (n probes = 64)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% MEST$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="PEG3 (n probes = 25)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% PEG3$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="NNAT (n probes = 9)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% NNAT$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="NDN (n probes = 6)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% NDN$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="SNRPN (n probes = 267)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% SNRPN$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="SGCE (n probes = 40)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% SGCE$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="PEG10 (n probes = 70)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% PEG10$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="MEG3 (n probes = 51)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% MEG3$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="H19 (n probes = 52)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% H19$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="PLAGL1 (n probes = 36)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% PLAGL1$name,])
imprinted.regions$ID[imprinted.regions$imprinted.region=="GRB10 (n probes = 52)"]<-1:nrow(imprinted.regions[imprinted.regions$MarkerName %in% GRB10$name,])

imprinted.regions$lower<- imprinted.regions$Effect-(1.96*imprinted.regions$StdErr)
imprinted.regions$upper<- imprinted.regions$Effect+(1.96*imprinted.regions$StdErr)
require(ggplot2)
require(wesanderson)

Plot<-ggplot(imprinted.regions, aes(x=ID))+
geom_ribbon(aes(ymin=lower*100,ymax=upper*100),alpha=0.3,fill=wes_palette("Zissou1")[1])+
geom_point(aes(y=Effect*100),size=0.5)+
geom_hline(yintercept=0)+
facet_wrap(~imprinted.region, ncol=3,scales="free")+
theme_grey(base_size=8) + labs(x = "Ordered position within region",y = "Difference in percentage methylation per 1 SD increase in paternal BMI",size=12) + 
theme(axis.text.y = element_text(colour="black",size=8),axis.text.x=element_blank())+
theme(strip.background = element_rect(fill="white"), strip.text=element_text(size=8), plot.background=element_rect(fill="white"))+
ggtitle("Associations between paternal BMI (adjusted for maternal BMI)\nand offspring methylation at birth in imprinted regions\n")+
theme(plot.title=element_text(size=14, hjust=0.5))


pdf("imprinted.patmatbmi.pdf")
Plot
dev.off()

## QQ plots

prep.myQQ <- function(p.sva, ci = 0.95) {
  n  <- length(p.sva)
  df <- data.frame(
    observed = -log10(sort(p.sva)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
}

X<-lapply(c("IGF2","MEST","PEG3","NNAT","NDN","SNRPN","SGCE","PEG10","MEG3","H19","PLAGL1","GRB10"),function(x) imprinted.regions[grep(imprinted.regions$imprinted.region,pattern=x),"Pvalue"])
X<- lapply(X,prep.myQQ)
X<-do.call(rbind,X)

X$imprinted.region <- c(rep("IGF2 (n probes = 97)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),])),
rep("MEST (n probes = 64)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),])),
rep("PEG3 (n probes = 25)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),])),
rep("NNAT (n probes = 9)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),])),
rep("NDN (n probes = 6)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),])),
rep("SNRPN (n probes = 267)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),])),
rep("SGCE (n probes = 40)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),])),
rep("PEG10 (n probes = 70)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),])),
rep("MEG3 (n probes = 51)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),])),
rep("H19 (n probes = 52)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),])),
rep("PLAGL1 (n probes = 36)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),])),
rep("GRB10 (n probes = 52)",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),]))
)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    facet_wrap(~imprinted.region, ncol=3,scales="free")+
    ggtitle("Paternal BMI (adjusted for maternal BMI) P-value distributions \nat imprinted regions") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  


pdf("imprinted.qq.patmatbmi.pdf")
Plot
dev.off()
          
## Kolmogorov-Smirnov Tests to see distribution of observed P-values fits to null distribution
          
X<-lapply(c("IGF2","MEST","PEG3","NNAT","NDN","SNRPN","SGCE","PEG10","MEG3","H19","PLAGL1","GRB10"),function(x) imprinted.regions[grep(imprinted.regions$imprinted.region,pattern=x),"Pvalue"])
ks.results<-unlist(lapply(lapply(X,function(x) ks.test(unlist(x),distribution="punif",0,1)),function(x) unlist(x)[2]))
  
                                 
#Donkin et al
## QQ at sites in the same direction
          
donkin <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/DONKIN.csv", header=TRUE)
donkin <- merge(donkin, annotation, by.x="start", by.y="position",all=F)
donkin <- merge(list.of.results$covs.patmat,donkin,by.x="MarkerName",by.y="name",all=F)
table(sign(donkin$Effect)==sign(donkin$meth.diff..Obese.versus.Lean.))
donkin <- donkin[sign(donkin$Effect)==sign(donkin$meth.diff..Obese.versus.Lean.),]
donkin$FDR <- p.adjust(donkin$Pvalue,method="fdr")

X<-prep.myQQ(donkin$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Paternal BMI (adjusted for maternal BMI) P-value distributions \nat regions identified in Donkin et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("donkin.qq.patmatbmi.pdf")
Plot
dev.off()
                                 
ks.test(donkin$Pvalue,distribution="punif",0,1)
                                 
#Sharp et al previous maternal BMI PACE study
## QQ patmat
          
sharp <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv", header=TRUE)
sharp <- merge(list.of.results$covs.patmat,sharp,by.x="MarkerName",by.y="CpG",all=F)#64
table(sign(sharp$Effect.x)==sign(sharp$Effect_cells))#28 same direction
sharp$FDR <- p.adjust(sharp$Pvalue,method="fdr")
summary(sharp$FDR<0.05)#0
                                 
X<-prep.myQQ(sharp$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Paternal BMI (adjusted for maternal BMI) P-value distributions \nat maternal BMI CpGs from Sharp et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("sharp.qq.patmatbmi.pdf")
Plot
dev.off()
                                 
ks.test(sharp$Pvalue,distribution="punif",0,1)#0.03077
                                 

## QQ pat
          
sharp <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv", header=TRUE)
sharp <- merge(list.of.results$covs.pat,sharp,by.x="MarkerName",by.y="CpG",all=F)#64
table(sign(sharp$Effect.x)==sign(sharp$Effect_cells))#43 same direction
sharp$FDR <- p.adjust(sharp$Pvalue,method="fdr")
summary(sharp$FDR<0.05)#0
                                 
X<-prep.myQQ(sharp$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Paternal BMI meta-EWAS P-value distributions \nat maternal BMI CpGs from Sharp et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("sharp.qq.patbmi.pdf")
Plot
dev.off()
                                 
ks.test(sharp$Pvalue,distribution="punif",0,1)#0.03077
                                 
                                 
 ## QQ patmat
          
sharp <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv", header=TRUE)
sharp <- merge(list.of.results$covs.patmat,sharp,by.x="MarkerName",by.y="CpG",all=F)#64
table(sign(sharp$Effect.x)==sign(sharp$Effect_cells))#28 same direction
sharp$FDR <- p.adjust(sharp$Pvalue,method="fdr")
summary(sharp$FDR<0.05)#0
                                 
X<-prep.myQQ(sharp$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Paternal BMI (adjusted for maternal BMI) meta-EWAS P-value distributions \nat maternal BMI CpGs from Sharp et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("sharp.qq.patmatbmi.pdf")
Plot
dev.off()
                                 
ks.test(sharp$Pvalue,distribution="punif",0,1)#0.03077
                                 

## QQ mat
          
sharp <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv", header=TRUE)
sharp <- merge(list.of.results$covs.mat,sharp,by.x="MarkerName",by.y="CpG",all=F)#64
table(sign(sharp$Effect.x)==sign(sharp$Effect_cells))#62 same direction
sharp$FDR <- p.adjust(sharp$Pvalue,method="fdr")
summary(sharp$FDR<0.05)#35
                                 
X<-prep.myQQ(sharp$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Maternal BMI meta-EWAS P-value distributions \nat maternal BMI CpGs from Sharp et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("sharp.qq.matbmi.pdf")
Plot
dev.off()
                                 
ks.test(sharp$Pvalue,distribution="punif",0,1)#0.03077
                                 
## QQ matpat
          
sharp <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv", header=TRUE)
sharp <- merge(list.of.results$covs.matpat,sharp,by.x="MarkerName",by.y="CpG",all=F)#64
table(sign(sharp$Effect.x)==sign(sharp$Effect_cells))#62 same direction
sharp$FDR <- p.adjust(sharp$Pvalue,method="fdr")
summary(sharp$FDR<0.05)#35
                                 
X<-prep.myQQ(sharp$Pvalue)

Plot <-  ggplot(X) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour=wes_palette("Zissou1")[5]) +
    geom_line(aes(expected, clower), linetype =2, colour=wes_palette("Zissou1")[5]) +
    geom_point(aes(expected, observed), shape = 20, size = 2,colour=wes_palette("Zissou1")[1]) +
    ggtitle("Maternal BMI (adjusted for paternal BMI) meta-EWAS P-value distributions \nat maternal BMI CpGs from Sharp et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("sharp.qq.matpatbmi.pdf")
Plot
dev.off()
                                 
ks.test(sharp$Pvalue,distribution="punif",0,1)#0.03077
