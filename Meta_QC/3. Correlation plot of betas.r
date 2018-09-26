# Cohort QC: Correlation plot of betas in all 12 results files

time_point <-"birth"#or whatever

require(corrplot)
require(plyr)

extract.coefficients <- function(ewas.dataframe){
	ewas.dataframe[,which(colnames(ewas.dataframe)=="Effect")]
}

correlation.plot<-function(cohort,cohort_name){
x <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))
colnames(x)<-names(list.of.results)
filename <- paste0(cohort_name,".correlation.",time_point,".png")
png(filename,width=15,height=18,units="cm",res=300)
corrplot(cor(x),method="number",type="upper")
title(cohort_name)
dev.off()
}

correlation.plot(list.of.results,"meta_models")
