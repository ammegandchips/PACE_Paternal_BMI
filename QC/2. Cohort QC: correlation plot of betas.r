# Cohort QC: Correlation plot of betas in all 12 results files

require(corrplot)
require(plyr)

extract.coefficients <- function(ewas.dataframe){
	ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")]
}

correlation.plot<-function(cohort,cohort_name){
x <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))
colnames(x)<-key$result[which(colnames(x)%in%key$merged.title.coef)]
filename <- paste0("qc_res/",cohort_name,".correlation.png")
png(filename,width=15,height=18,units="cm",res=300)
corrplot(cor(x),method="number",type="upper")
title(cohort_name)
dev.off()
}

lapply(1:length(All.EWAS),function(x) correlation.plot(cohort=All.EWAS[[x]],cohort_name=names(All.EWAS)[[x]]))
