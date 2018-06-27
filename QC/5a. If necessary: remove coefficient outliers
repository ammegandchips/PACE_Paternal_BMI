#Remove outliers in effect sizes (anything with e.s. outside the 99.999%ile)

#Finding cut-off:

require(reshape)

abbreviated.cohort.names <- c("AL","Bas","Bwh","CH","EN","GR","Inc","Ic","PI","VI","RH")

list.of.dataframes <- lapply(All.EWAS,
function(cohort){
x <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))
colnames(x)<-key$result[which(colnames(x)%in%key$merged.title.coef)]
x
})

list.of.melted.dataframes <- lapply(list.of.dataframes,melt)
names(list.of.melted.dataframes) <- abbreviated.cohort.names
dat <- data.frame(do.call(rbind, list.of.melted.dataframes))
dat$Cohorts <- row.names(dat)
dat$Cohorts <- unlist(lapply(strsplit(dat$Cohorts,split=".",fixed=TRUE),"[",1))
colnames(dat) <- c("Models","Coefficients","Cohorts")
cutoff.lower <- quantile(dat$Coefficients,1e-5,na.rm=T) 
cuttoff.higher <- quantile(dat$Coefficients,1-1e-5,na.rm=T)
cutoff <- round(max(abs(cutoff.lower),abs(cutoff.higher)),1)

remove.outlying.es <- function(ewas.dataframe){
  cbind(ewas.dataframe[,-grep(colnames(ewas.dataframe),pattern="coef")],
    replace(ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")],abs(ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")])>cutoff,NA))
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,remove.outlying.es))

#The above function gives the affected columns a weird name, so this function renames them with something more sensible:
replace.name<-function(X){
  colnames(X)[grep(colnames(X),pattern="ewas")]<-"coef"
  X
}

All.EWAS <- lapply(All.EWAS,function(x) lapply(x,replace.name))
