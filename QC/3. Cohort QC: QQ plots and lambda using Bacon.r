# Cohort QC: QQ plots and lambda using Bacon

time_point <- "birth" #or whatever

require(bacon)
require(gridExtra)
extract.se <- function(ewas.dataframe){
	ewas.dataframe[,c(which(colnames(ewas.dataframe)=="se"),grep(colnames(ewas.dataframe),pattern="se.pheno"))]
}

bacon.plot <- function(cohort,cohort_name){
z.scores <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))/data.frame(do.call(cbind, lapply(cohort,extract.se)))
colnames(z.scores)<-key$result[which(colnames(z.scores)%in%key$merged.title.coef)]
bacon.res <- apply(z.scores,2,bacon)
bacon.lambdas <- data.frame(do.call(cbind,lapply(bacon.res,inflation)))
filename <- paste0("qc_res/",cohort_name,".",time_point,".qqplots.png")
png(filename, width=30,height=40,units="cm",res=300)
plot.function<-function(i){
	plot(bacon.res[[i]],type="qq")+ggtitle(paste(cohort_name,names(bacon.res)[[i]],"lambda =",round(bacon.lambdas[[i]],3)))+theme(legend.position="none",plot.title=element_text(hjust=0.5),axis.title=element_text(size=8))
}
plots <- lapply(1:length(bacon.res),plot.function)
do.call(grid.arrange,plots)
dev.off()
}

lapply(1:length(All.EWAS),function(x) bacon.plot(cohort=All.EWAS[[x]],cohort_name=names(All.EWAS)[[x]]))
