#Volcano plots

require(ggplot2)

Pat.volcano<-ggplot(list.of.results$covs.pat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Paternal BMI meta-EWAS results at birth")

Mat.volcano<-ggplot(list.of.results$covs.mat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Maternal BMI meta-EWAS results at birth")

PatMat.volcano<-ggplot(list.of.results$covs.patmat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Paternal BMI (adjusted for maternal BMI) \nmeta-EWAS results at birth")

MatPat.volcano<-ggplot(list.of.results$covs.matpat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Maternal BMI (adjusted for paternal BMI) \nmeta-EWAS results at birth")

PatMat.male.volcano<-ggplot(list.of.results$boys.patmat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Paternal BMI adjusted for maternal BMI \n(male offspring) meta-EWAS results at birth")

MatPat.male.volcano<-ggplot(list.of.results$boys.matpat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+

ggtitle("Maternal BMI adjusted for paternal BMI \nmeta-EWAS results at birth (male offspring)")

PatMat.female.volcano<-ggplot(list.of.results$girls.patmat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+

ggtitle("Paternal BMI adjusted for maternal BMI \nmeta-EWAS results at birth (female offspring) ")

MatPat.female.volcano<-ggplot(list.of.results$girls.matpat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5,size=Pvalue<1e-5))+
scale_colour_manual(values=c("grey14","#3B9AB2"))+
scale_size_manual(values=c(1,2))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7.5))+
xlim(c(-1.1,1.1))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5,size=16),axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
ggtitle("Maternal BMI adjusted for paternal BMI \nmeta-EWAS results at birth (female offspring)")


png("Volcano.patbmi.png",width=500,height=500)
Pat.volcano
dev.off()

png("Volcano.matbmi.png",width=500,height=500)
Mat.volcano
dev.off()

png("Volcano.patmatbmi.png",width=500,height=500)
PatMat.volcano
dev.off()

png("Volcano.matpatbmi.png",width=500,height=500)
MatPat.volcano
dev.off()

png("Volcano.patmatbmi.male.png",width=500,height=500)
PatMat.male.volcano
dev.off()

png("Volcano.matpatbmi.male.png",width=500,height=500)
MatPat.male.volcano
dev.off()

png("Volcano.patmatbmi.female.png",width=500,height=500)
PatMat.female.volcano
dev.off()

png("Volcano.matpatbmi.female.png",width=500,height=500)
MatPat.female.volcano
dev.off()
