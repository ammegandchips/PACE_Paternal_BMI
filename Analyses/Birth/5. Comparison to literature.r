#In the literature, we found evidence for association with paternal BMI at imprinted regions and at >9000 CpGs in a paper by Donkin et al.

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

list.of.results$covs.patmat <- list.of.results$covs.patmat[order(list.of.results$covs.patmat$Effect),]

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

imprinted.regions$ID <- c(1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),]))

require(ggplot2)

Plot<-ggplot(imprinted.regions, aes(x=1, y=ID))+
geom_tile(aes(fill = Effect))+
facet_wrap(~imprinted.region, ncol=3,scales="free")+
xlab("")+ylab("")+
scale_fill_gradientn(colours = c("darkblue", "white", "red"),name="Effect\nestimate",limits=c(-max(abs(imprinted.regions$Effect)),max(abs(imprinted.regions$Effect))))+
theme_grey(base_size=8) + labs(x = "",y = "") + 
scale_x_discrete(expand = c(0, 0)) +
scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(),axis.text.y = element_text(colour="black",size=11))+
theme(strip.background = element_rect(fill="grey90"), strip.text=element_text(size=12), plot.background=element_rect(fill="grey90"))+
ggtitle("Associations between paternal BMI and \noffspring cord blood methylation at imprinted regions\n")+
theme(plot.title=element_text(size=16, hjust=0.5))+
theme(legend.position="bottom",legend.text = element_text(size = 12),legend.title=element_text(size=12),legend.background=(element_rect(fill=NA)),legend.key.width=unit(2.5,"cm"))


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
    geom_point(aes(expected, observed), shape = 20, size = 2,colour="cornflowerblue") +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    facet_wrap(~imprinted.region, ncol=3,scales="free")+
    geom_line(aes(expected, clower), linetype = 2) +
    ggtitle("P-value enrichment at imprinted regions") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="grey90"),panel.background=element_rect(fill="grey90"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  


pdf("imprinted.qq.patmatbmi.pdf")
Plot
dev.off()
          
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
    geom_point(aes(expected, observed), shape = 20, size = 2,colour="cornflowerblue") +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    ggtitle("P-value enrichment at regions identified in Donkin et al.") +
  theme_classic()+
    theme(plot.title=element_text(size=16, hjust=0.5),strip.background=element_blank(),strip.text=element_text(size=12),plot.background=element_rect(fill="grey90"),panel.background=element_rect(fill="grey90"))+
    xlab(expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))  

pdf("donkin.qq.patmatbmi.pdf")
Plot
dev.off()
