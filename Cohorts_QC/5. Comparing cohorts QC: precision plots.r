# Comparing cohorts QC: Precision plot (1/median[SE] vs sqrt(N)) for each model

require(data.table)

abbreviated.cohort.names <- c("AL","Bas","Bwh","CH","EN","GR","GO","Inc","Ic","M1","M2","M3","PI","PV","RH")
#abbreviated.cohort.names <- c("AL","CH","GR","Inc","Ic","PV") #for childhood
time_point <- "birth" #or childhood

extract.median.n <- function(ewas.dataframe){
	median(ewas.dataframe[,c(which(colnames(ewas.dataframe)=="n"))],na.rm=TRUE)
}

list.of.vectors <- lapply(All.EWAS,
function(cohort){
do.call(cbind, lapply(cohort,extract.median.n))
})

names(list.of.vectors) <- abbreviated.cohort.names
sample.sizes <- as.data.frame(transpose(list.of.vectors),row.names=names(list.of.vectors),col.names=colnames(list.of.vectors[[1]]))
colnames(sample.sizes) <- c("min.pat","min.mat","min.patmat","covs.pat","covs.mat","covs.patmat","boys.patmat","girls.patmat")
sample.sizes$min.matpat <- sample.sizes$min.patmat
sample.sizes$covs.matpat <- sample.sizes$covs.patmat
sample.sizes$girls.matpat <- sample.sizes$girls.patmat
sample.sizes$boys.matpat <- sample.sizes$boys.patmat
sample.sizes.melted <- melt(sample.sizes)
sample.sizes.melted$cohort<-row.names(sample.sizes)
names(sample.sizes.melted) <- c("model","sample.size","cohort")

# head(sample.sizes.melted)
#    model sample.size cohort
# min.pat         647     AL
# min.pat          88    Bas
# min.pat         138    Bwh
# min.pat         170     CH
# min.pat         183     EN
# min.pat        1061     GR

extract.se <- function(ewas.dataframe){
	ewas.dataframe[,c(which(colnames(ewas.dataframe)=="se"),grep(colnames(ewas.dataframe),pattern="se.pheno"))]
}

list.of.dataframes <- lapply(All.EWAS,
function(cohort){
x <- data.frame(do.call(cbind, lapply(cohort,extract.se)))
colnames(x)<-key$result[which(colnames(x)%in%key$merged.title.se)]
x
})

list.of.medians <- lapply(list.of.dataframes,function(x) apply(x,2,median,na.rm=TRUE))
names(list.of.medians) <- abbreviated.cohort.names
median.SEs <- as.data.frame(transpose(list.of.medians),row.names=names(list.of.medians),col.names=names(list.of.medians[[1]]))
median.SEs.melted <- melt(median.SEs)
median.SEs.melted$cohort<-row.names(median.SEs)
names(median.SEs.melted) <- c("model","median.se","cohort")

# head(median.SEs.melted)
#    model    median.se cohort
# min.pat 0.0010936286     AL
# min.pat 0.0025318482    Bas
# min.pat 0.0022816048    Bwh
# min.pat 0.0010440241     CH
# min.pat 0.0008717286     EN
# min.pat 0.0005681995     GR

dat <- cbind(sample.sizes.melted,median.SEs.melted$median.se)
names(dat) <- c("Model","N","Cohort","Median.SE")

#    Model    N Cohort    Median.SE
# min.pat  647     AL 0.0010936286
# min.pat   88    Bas 0.0025318482
# min.pat  138    Bwh 0.0022816048
# min.pat  170     CH 0.0010440241
# min.pat  183     EN 0.0008717286
# min.pat 1061     GR 0.0005681995

require(ggplot2)
filename <- paste0("qc_res/Precision.by.samplesize.",time_point,".png")
png(filename, width=30,height=40,units="cm",res=300)
P <- ggplot(dat, aes(x = sqrt(N), y = 1/Median.SE, label=Cohort)) +
        geom_text()+
        ylab("Precision (1/median(SE))")+
        xlab("sqrt(sample size)") +
        ggtitle("Precision by Sample Size") +
        theme_bw() +
        theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold", hjust=0.5),
              text = element_text(size = 12, family = "Tahoma"),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 11)) +
        facet_wrap( ~ Model,ncol=3,scales="free")
P
dev.off()
