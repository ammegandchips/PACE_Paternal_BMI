time_point <- "birth" #or childhood

require(gridExtra)

myQQ <- function(P, ci = 0.95,Title) {
  lambda <- Lambda(P)
  n  <- length(P)
  df <- data.frame(
    observed = -log10(P),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,colour="blue") +
#	ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2,colour="blue") +
    ggtitle(paste0(Title," Lambda = ",round(lambda,2))) +
	theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(log10Pe) +
    ylab(log10Po)  
} 

Lambda<-function(P){
chisq <- qchisq(1-P,1)
median(chisq,na.rm=T)/qchisq(0.5,1)
}

extract.p <- function(ewas.dataframe){
	ewas.dataframe[,c(which(colnames(ewas.dataframe)=="p"),grep(colnames(ewas.dataframe),pattern="p.pheno"))]
}

qq.plot <- function(cohort,cohort_name){
Ps <- data.frame(do.call(cbind, lapply(cohort,extract.p)))
colnames(Ps)<-key$result[which(colnames(Ps)%in%key$merged.title.p)]
filename <- paste0("qc_res/",cohort_name,".",time_point,".qqplots.png")
png(filename, width=30,height=40,units="cm",res=300)
plot.function<-function(i){
	myQQ(sort(Ps[,i]),Title=paste(cohort_name,": ",colnames(Ps)[i]))
	}
plots <- lapply(1:ncol(Ps),plot.function)
do.call(grid.arrange,plots)
dev.off()
}

lapply(1:length(All.EWAS),function(x) qq.plot(cohort=All.EWAS[[x]],cohort_name=names(All.EWAS)[[x]]))
