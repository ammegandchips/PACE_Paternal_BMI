#agreement and correlation of pat vs patmat
cor(pat.man$Effect,patmat.man$Effect)#0.97  

plot(pat.man$Effect,patmat.man$Effect,main="Correlation between paternal BMI effect estimates \nobtained before and after adjustment for maternal BMI: r=0.97",xlab="Unadjusted",ylab="Adjusted")
abline(0,1,col=wes_palette("Zissou1")[1])

summary(abs(abs(patmat.man$Effect*100 - pat.man$Effect*100)/(pat.man$Effect*100)))

#comparing pat vs mat effect sizes
summary(abs(list.of.results$covs.patmat$Effect)< abs(list.of.results$covs.matpat$Effect))
summary(abs(list.of.results$covs.pat$Effect)< abs(list.of.results$covs.mat$Effect))

#Pat vs Mat manhattan plots

require(meffil)
require(wesanderson)

desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

annotation <- meffil.get.features("450k")
source("~/EWAS/pat_bmi/manhattan_plot_function.r")

pat.man <- merge(list.of.results$covs.pat,annotation,by.x="MarkerName",by.y="name",all=F)
mat.man <- merge(list.of.results$covs.mat,annotation,by.x="MarkerName",by.y="name",all=F)
pat.man$CHR <- as.numeric(unlist(lapply(pat.man$chromosome,substring,first=4)))
mat.man$CHR <- as.numeric(unlist(lapply(mat.man$chromosome,substring,first=4)))
patmat.man <- merge(list.of.results$covs.patmat,annotation,by.x="MarkerName",by.y="name",all=F)
matpat.man <- merge(list.of.results$covs.matpat,annotation,by.x="MarkerName",by.y="name",all=F)
patmat.man$CHR <- as.numeric(unlist(lapply(patmat.man$chromosome,substring,first=4)))
matpat.man$CHR <- as.numeric(unlist(lapply(matpat.man$chromosome,substring,first=4)))

pat.man$abs.Effect <-abs(pat.man$Effect)*100
mat.man$abs.Effect <-abs(mat.man$Effect)*100
patmat.man$abs.Effect <-abs(patmat.man$Effect)*100
matpat.man$abs.Effect <-abs(matpat.man$Effect)*100

png("pat.manhattan.effect.png",width=1500,height=500)
manhattan (pat.man, chr = "CHR", bp = "position", p = "abs.Effect", snp = "MarkerName", col = c(wes_palette("Zissou1")[1], 
    desat(wes_palette("Zissou1")[1])), chrlabs = NULL, suggestiveline = NULL, 
    genomewideline = NULL,ylim_max=max(c(pat.man$abs.Effect,mat.man$abs.Effect,matpat.man$abs.Effect,patmat.man$abs.Effect))+0.2,logp=F)
dev.off()

png("mat.manhattan.effect.png",width=1500,height=500)
manhattan (mat.man, chr = "CHR", bp = "position", p = "abs.Effect", snp = "MarkerName", col = c(wes_palette("Zissou1")[5], 
    desat(wes_palette("Zissou1")[5])), chrlabs = NULL, suggestiveline = NULL, 
    genomewideline = NULL,ylim_max=max(c(pat.man$abs.Effect,mat.man$abs.Effect,matpat.man$abs.Effect,patmat.man$abs.Effect))+0.2,logp=F)
dev.off()

png("patmat.manhattan.effect.png",width=1500,height=500)
manhattan (patmat.man, chr = "CHR", bp = "position", p = "abs.Effect", snp = "MarkerName", col = c(wes_palette("Zissou1")[1], 
    desat(wes_palette("Zissou1")[1])), chrlabs = NULL, suggestiveline = NULL, 
    genomewideline = NULL,ylim_max=max(c(pat.man$abs.Effect,mat.man$abs.Effect,matpat.man$abs.Effect,patmat.man$abs.Effect))+0.2,logp=F)
dev.off()

png("matpat.manhattan.effect.png",width=1500,height=500)
manhattan (matpat.man, chr = "CHR", bp = "position", p = "abs.Effect", snp = "MarkerName", col = c(wes_palette("Zissou1")[5], 
    desat(wes_palette("Zissou1")[5])), chrlabs = NULL, suggestiveline = NULL, 
    genomewideline = NULL,ylim_max=max(c(pat.man$abs.Effect,mat.man$abs.Effect,matpat.man$abs.Effect,patmat.man$abs.Effect))+0.2,logp=F)
dev.off()

#mat vs pat meta-analysis

fixed.effects.meta.analysis <- function(data){
                              require(metafor)
                              res = split(data, f=data$MarkerName)
                              res = lapply(res, function(x) rma.uni(slab=x$study,yi=x$Effect,sei=x$StdErr,method="FE",weighted=TRUE))
                              res
                              }

CpGs <- list.of.results$covs.pat[which(list.of.results$covs.pat$Pvalue<1e-5),"MarkerName"]

PatMatComparison<-rbind(list.of.results$covs.pat[which(list.of.results$covs.pat$MarkerName %in% CpGs),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% CpGs),],
list.of.results$covs.mat[which(list.of.results$covs.mat$MarkerName %in% CpGs),],
list.of.results$covs.matpat[which(list.of.results$covs.matpat$MarkerName %in% CpGs),])
PatMatComparison$Model <- c(rep("Paternal",length(CpGs)),rep("Paternal adjusted for maternal",length(CpGs)),
						rep("Maternal",length(CpGs)),rep("Maternal adjusted for paternal",length(CpGs)))
							  
results.patmat.adj<-fixed.effects.meta.analysis(data=PatMatComparison[which(PatMatComparison$Model %in% c("Paternal adjusted for maternal","Maternal adjusted for paternal")),])
results.patmat<-ldply(lapply(results.patmat,function(x) unlist(c(x[c("QE","QEp","I2")]))))
results.patmat<-fixed.effects.meta.analysis(data=PatMatComparison[which(PatMatComparison$Model %in% c("Paternal","Maternal")),])
results.patmat.adj<-ldply(lapply(results.patmat.adj,function(x) unlist(c(x[c("QE","QEp","I2")]))))
write.csv(results.patmat.adj,"matpatvspatmat.metaanalysis.birth.csv")

PatMatComparison$colour <-wes_palette("Zissou1")[1]
PatMatComparison$colour[which(PatMatComparison$Model %in% c("Maternal","Maternal adjusted for paternal"))]<-wes_palette("Zissou1")[5]
PatMatComparison$ci.lb<-PatMatComparison$Effect - (1.96* PatMatComparison$StdErr)
PatMatComparison$ci.ub<-PatMatComparison$Effect + (1.96* PatMatComparison$StdErr)
PatMatComparison$Model<-factor(PatMatComparison$Model,levels=c("Paternal","Paternal adjusted for maternal","Maternal","Maternal adjusted for paternal"),ordered=TRUE)
Order<-list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% CpGs),]
Order<-Order[order(Order$Effect),]
PatMatComparison$MarkerName<-factor(PatMatComparison$MarkerName,levels=Order$MarkerName,ordered=TRUE)
PatMatComparison<-merge(PatMatComparison,annotation,by.x="MarkerName",by.y="name",all=FALSE)
new.annotations<-fread("/panfs/panasas01/sscm/gs8094/Common_files/enhanced_annotations.txt",stringsAsFactors=FALSE)
PatMatComparison<-merge(PatMatComparison,new.annotations,by.x="MarkerName",by.y="IlmnID",all.y=FALSE)
PatMatComparison$gene <- unlist(lapply(strsplit(PatMatComparison$gene.symbol,split=";"),"[",1))
PatMatComparison$gene[is.na(PatMatComparison$gene)]<-PatMatComparison$"UCSC KnownGene"[is.na(PatMatComparison$gene)]
PatMatComparison$CpG.Gene <- paste0(PatMatComparison$MarkerName,"\n",PatMatComparison$gene)
PatMatComparison<-PatMatComparison[order(PatMatComparison$MarkerName),]
PatMatComparison$CpG.Gene<-factor(PatMatComparison$CpG.Gene,levels=unique(PatMatComparison$CpG.Gene),ordered=TRUE)

P <- ggplot(PatMatComparison,aes(x=Model,y=Effect*100))+
geom_hline(yintercept=0,linetype="dashed")+
geom_errorbar(aes(colour=Model,ymin=ci.lb*100, ymax=ci.ub*100),width=0.5,size=1)+
geom_point(aes(shape=Model,colour=Model),fill="white",size=4)+
scale_shape_manual(values=c(15,22,19,21))+
scale_colour_manual(values=c(wes_palette("Zissou1")[1],wes_palette("Zissou1")[1],wes_palette("Zissou1")[5],wes_palette("Zissou1")[5]))+
facet_grid(.~CpG.Gene)+
theme_bw() + theme(legend.spacing.x = unit(1.0, 'cm'),legend.position = "bottom",legend.title=element_blank(),legend.text=element_text(size=14),axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14),panel.grid.major.x = element_blank()) +
xlab("")+ylab("Effect estimate (difference in % methylation\nper 1SD increase in parental BMI)")+
ggtitle("CpGs showing strongest statistical evidence of association with paternal BMI\n") + theme(panel.spacing = unit(0.8, "lines"),panel.border = element_blank(),panel.background = element_rect(fill="grey95"),plot.title = element_text(hjust = 0.5,size=16),strip.background=element_rect(fill="grey95",colour=NA),strip.text = element_text(size=12,face = "italic"))

				 
png("PatvsMat.coefplot.png",width=1000,height=500)
P
dev.off()
				 
