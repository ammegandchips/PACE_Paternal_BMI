#Leave one out analysis for top sites (P<1e-5) from covs adjusted models

#E.g. using paternal models

top.sites.covs.pat <- list.of.results$covs.pat[which(list.of.results$covs.pat$Pvalue<1e-5),"MarkerName"]
top.sites.covs.patmat <- list.of.results$covs.patmat[which(list.of.results$covs.patmat$Pvalue<1e-5),"MarkerName"]
top.sites.boys <- list.of.results$boys.patmat[which(list.of.results$boys.patmat$Pvalue<1e-5),"MarkerName"]
top.sites.girls <- list.of.results$girls.patmat[which(list.of.results$girls.patmat$Pvalue<1e-5),"MarkerName"]

addprobeID <- function(ewas.dataframe){
  ewas.dataframe$probeid <- row.names(ewas.dataframe)
  ewas.dataframe
}

All.EWAS.top <- lapply(All.EWAS,function(x) lapply(x,addprobeID))

All.EWAS.top.covs.pat<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.covs.pat),]))
All.EWAS.top.covs.patmat<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.covs.patmat),]))
All.EWAS.top.boys<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.boys),]))
All.EWAS.top.girls<-lapply(All.EWAS.top,function(x) lapply(x,function(y) y[which(y$probeid %in% top.sites.girls),]))

covs.pat<-lapply(All.EWAS.top.covs.pat,function(x) x$ewas.res.covs.mutual)
covs.patmat<-lapply(All.EWAS.top.covs.patmat,function(x) x$ewas.res.covs.mutual)
boys<-lapply(All.EWAS.top.boys,function(x) x$ewas.res.covs.mutual.boys.only)
girls<-lapply(All.EWAS.top.girls,function(x) x$ewas.res.covs.mutual.girls.only)

# Remove cohorts if they didn't contribute results

covs.pat$ENVIRONAGE <- NULL
covs.patmat$ENVIRONAGE <- NULL
boys$ENVIRONAGE <- NULL
girls$ENVIRONAGE <- NULL
girls$BIB_asian<- NULL

# Prepare for meta-analysis (necessary for L-O-O)

covs.pat<-do.call(rbind,covs.pat)
covs.patmat<-do.call(rbind,covs.patmat)
boys<-do.call(rbind,boys)
girls<-do.call(rbind,girls)

covs.pat$study <- unlist(lapply(strsplit(row.names(covs.pat),split=".",fixed=TRUE),"[",1))
covs.patmat$study <- unlist(lapply(strsplit(row.names(covs.patmat),split=".",fixed=TRUE),"[",1))
boys$study <- unlist(lapply(strsplit(row.names(boys),split=".",fixed=TRUE),"[",1))
girls$study <- unlist(lapply(strsplit(row.names(girls),split=".",fixed=TRUE),"[",1))

fixed.effects.meta.analysis <- function(data,coef.name,se.name){
                              require(metafor)
                              data <- data[,c("probeid","study",coef.name,se.name)]
                              colnames(data) <- c("probeid","study","coef","se")
                              res = split(data, f=data$probeid)
                              res = lapply(res, function(x) rma.uni(slab=x$study,yi=x$coef,sei=x$se,method="FE",weighted=TRUE))
                              res
                              }

#Run meta-analysis

covs.pat<-fixed.effects.meta.analysis(data=covs.pat,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
covs.patmat<-fixed.effects.meta.analysis(data=covs.patmat,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
boys<-fixed.effects.meta.analysis(data=boys,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")
girls<-fixed.effects.meta.analysis(data=girls,coef.name="coef.pheno.data.Zpat.bmi",se.name="se.pheno.data.Zpat.bmi")

#Run L-O-O

covs.pat.loo <- lapply(lapply(covs.pat,leave1out),function(x) as.data.frame(print(x)))
covs.patmat.loo <- lapply(lapply(covs.patmat,leave1out),function(x) as.data.frame(print(x)))
boys.loo <- lapply(lapply(boys,leave1out),function(x) as.data.frame(print(x)))
girls.loo <- lapply(lapply(girls,leave1out),function(x) as.data.frame(print(x)))

#Combine original (full) results with L-O-O results

covs.pat.original = lapply(covs.pat, function(x)unlist(c(x$b[[1]],x[c("se","zval","pval","ci.lb","ci.ub","QE","QEp")])))
covs.patmat.original = lapply(covs.patmat, function(x)unlist(c(x$b[[1]],x[c("se","zval","pval","ci.lb","ci.ub","QE","QEp")])))
girls.original = lapply(girls, function(x)unlist(c(x$b[[1]],x[c("se","zval","pval","ci.lb","ci.ub","QE","QEp")])))
boys.original = lapply(boys, function(x)unlist(c(x$b[[1]],x[c("se","zval","pval","ci.lb","ci.ub","QE","QEp")])))

loo.processing<-
function(i,X.loo,X.original){
  res = rbind(X.loo[[i]],X.original[[i]])
  res$leftout<-factor(c(row.names(X.loo[[i]]),"none"),levels=c("none",sort(row.names(X.loo[[i]])),ordered=TRUE))
  res$colour <- c(rep("black",nrow(X.loo[[i]])),"red")
  res
}

covs.pat.loo <-lapply(1:length(covs.pat.loo),loo.processing,X.loo=covs.pat.loo,X.original=covs.pat.original)
covs.patmat.loo <-lapply(1:length(covs.patmat.loo),loo.processing,X.loo=covs.patmat.loo,X.original=covs.patmat.original)
boys.loo <-lapply(1:length(boys.loo),loo.processing,X.loo=boys.loo,X.original=boys.original)
girls.loo <-lapply(1:length(girls.loo),loo.processing,X.loo=girls.loo,X.original=girls.original)

names(covs.pat.loo) <- names(covs.pat)
names(covs.patmat.loo) <- names(covs.patmat)
names(boys.loo) <- names(boys)
names(girls.loo) <- names(girls)

# Draw coefficient plots for L-O-O analysis

require(ggplot2)

draw.loo.plot<- function(list.of.dfs,i){
  df <- list.of.dfs[[i]]
  name.df <- names(list.of.dfs)[[i]]
  ggplot(df,aes(x=leftout,y=estimate,colour=colour))+
geom_point()+
geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub),width=0.5)+
scale_color_manual(values=c("red","black"))+
geom_hline(yintercept=0,linetype="dashed")+
theme_bw() + theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1),panel.grid.major.x = element_blank()) +
xlab("Left out study")+ylab("Effect estimate (difference in methylation\nper 1SD increase in paternal BMI)")+
ggtitle(name.df) + theme(plot.title = element_text(hjust = 0.5))
}

pdf("looplots.covs.pat.birth.pdf",width=6,height=4)
lapply(1:length(covs.pat.loo),draw.loo.plot,list.of.dfs=covs.pat.loo)
dev.off()

pdf("looplots.covs.patmat.birth.pdf",width=6,height=4)
lapply(1:length(covs.patmat.loo),draw.loo.plot,list.of.dfs=covs.patmat.loo)
dev.off()

pdf("looplots.boys.birth.pdf",width=6,height=4)
lapply(1:length(boys.loo),draw.loo.plot,list.of.dfs=boys.loo)
dev.off()

pdf("looplots.girls.birth.pdf",width=6,height=4)
lapply(1:length(girls.loo),draw.loo.plot,list.of.dfs=girls.loo)
dev.off()
