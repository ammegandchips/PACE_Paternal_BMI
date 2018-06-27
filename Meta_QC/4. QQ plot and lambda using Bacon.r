#QQ plot and Lambda using Bacon for each model

require(bacon)
require(gridExtra)

extract.se <- function(ewas.dataframe){
	ewas.dataframe[,which(colnames(ewas.dataframe)=="StdErr")]
}

bacon.plot <- function(cohort,cohort_name){
z.scores <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))/data.frame(do.call(cbind, lapply(cohort,extract.se)))
bacon.res <- apply(z.scores,2,bacon)
bacon.lambdas <- data.frame(do.call(cbind,lapply(bacon.res,inflation)))
filename <- paste0(cohort_name,".qqplots.outliersremoved.png")
png(filename, width=30,height=40,units="cm",res=300)
plot.function<-function(i){
	plot(bacon.res[[i]],type="qq")+ggtitle(paste(cohort_name,names(bacon.res)[[i]],"lambda =",round(bacon.lambdas[[i]],3)))+theme(legend.position="none",plot.title=element_text(hjust=0.5),axis.title=element_text(size=8))
}
plots <- lapply(1:length(bacon.res),plot.function)
do.call(grid.arrange,plots)
dev.off()
}

bacon.plot(list.of.results,"meta_models")
