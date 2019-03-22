# Comparing cohorts QC: Boxplot to compare the distribution of betas between cohorts for each of the 12 models 
# Check no outliers, all in correct range, etc

require(reshape)

extract.coefficients <- function(ewas.dataframe){
	ewas.dataframe[,grep(colnames(ewas.dataframe),pattern="coef")]
}

abbreviated.cohort.names <- c("AL","Bas","Bwh","CH","GR","GO","IN","Mo1","Mo2","Mo3","PI","PV","RH")#or c("AL","CH","GR","HE","IN", PV") for childhood
time_point <- "birth" #or whatever

list.of.dataframes <- lapply(All.EWAS,
function(cohort){
x <- data.frame(do.call(cbind, lapply(cohort,extract.coefficients)))
colnames(x)<-key$result[which(colnames(x)%in%key$merged.title.coef)]
x
}
)

list.of.melted.dataframes <- lapply(list.of.dataframes,melt)
names(list.of.melted.dataframes) <- abbreviated.cohort.names
dat <- data.frame(do.call(rbind, list.of.melted.dataframes))
dat$Cohorts <- row.names(dat)
dat$Cohorts <- unlist(lapply(strsplit(dat$Cohorts,split=".",fixed=TRUE),"[",1))
colnames(dat) <- c("Models","Coefficients","Cohorts")

require(ggplot2)

filename <- paste0("qc_res/coefficient.distributions.",time_point,".png")
png(filename, width=30,height=40,units="cm",res=300)
P <- ggplot(dat, aes(x = Cohorts, y = Coefficients)) +
        geom_boxplot(fill = "#4271AE", colour = "#1F3552",
                     alpha = 0.7) +
        ylab("Effect estimate\n(difference in methylation per 1SD increase in parental BMI)")+
        xlab("Cohort") +
        ggtitle("Distribution of effect estimates") +
        theme_bw() +
        theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold", hjust=0.5),
              text = element_text(size = 12, family = "Tahoma"),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 11)) +
        facet_wrap( ~ Models,ncol=3,scales="free")
P
dev.off()
