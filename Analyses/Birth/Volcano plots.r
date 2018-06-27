#Volcano plots

require(ggplot2)

Pat.volcano<-ggplot(list.of.results.het.removed$covs.pat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5),size=0.4)+
scale_colour_manual(values=c("grey22","deepskyblue2"))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7))+
xlim(c(-4,4))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
ggtitle("Paternal BMI")

Mat.volcano<-ggplot(list.of.results.het.removed$covs.mat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5),size=0.4)+
scale_colour_manual(values=c("grey22","deepskyblue2"))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7))+
xlim(c(-4,4))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
ggtitle("Maternal BMI")

PatMat.volcano<-ggplot(list.of.results.het.removed$covs.patmat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5),size=0.4)+
scale_colour_manual(values=c("grey22","deepskyblue2"))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in paternal BMI")+
ylim(c(0,7))+
xlim(c(-4,4))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
ggtitle("Paternal BMI adjusted for maternal BMI")

MatPat.volcano<-ggplot(list.of.results.het.removed$covs.matpat,aes(x=Effect*100,y=-log10(Pvalue)))+
geom_point(aes(col=Pvalue<1e-5),size=0.4)+
scale_colour_manual(values=c("grey22","deepskyblue2"))+
theme_classic()+
xlab("Difference in % offspring DNA methylation\n per 1SD increase in maternal BMI")+
ylim(c(0,7))+
xlim(c(-4,4))+
ylab("-log10(P-value)")+
theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
ggtitle("Maternal BMI adjusted for paternal BMI")

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
