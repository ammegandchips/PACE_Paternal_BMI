#Read results for cell types ~ paternal BMI
ALSPAC <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/alspac/results/ALSPAC.patbmi.cells.res.summary.late_childhood.csv",stringsAsFactors=FALSE)
GenerationR <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/GenR/GenR.patbmi.cells.res.summary.late_childhood.csv",stringsAsFactors=FALSE)
CHAMACOS <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CHAMACOS_results/CHAMACOS.patbmi.cells.res.summary.late_childhood.csv",stringsAsFactors=FALSE)
HELIX <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/HELIX/HELIX.patbmi.cells.res.summary.late_childhood.csv",stringsAsFactors=FALSE)
INMA <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/results_PACE_INMA/4years_nocombat/PACE.patbmi.cells.res.summary.early_childhood.csv",stringsAsFactors=FALSE)
ProjectViva <- read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/VIVA_PATBMI/VIVA_PATBMI/Viva.patbmi.cells.res.summary.early_childhood.csv",stringsAsFactors=FALSE)


cell.results<-list(ALSPAC,CHAMACOS,GenerationR,HELIX,INMA,ProjectViva)
names(cell.results) <-c("ALSPAC","CHAMACOS","GenerationR","HELIX","INMA","ProjectViva")

cell.results <- lapply(cell.results,setNames, c("cell_type","effect","se","t","p"))
cell.results <- do.call(cbind,cell.results)

#meta-analysis

require(metafor)

studies <-c("ALSPAC","CHAMACOS","GenerationR","HELIX","INMA","ProjectViva")

fixed.effects.meta.analysis <- function(list.of.studies,data){
                              coefs = data[,c("ALSPAC.cell_type",paste0(list.of.studies,".effect"))]
							  colnames(coefs) <- c("cell",paste0(list.of.studies))
                              ses = data[,c("ALSPAC.cell_type",paste0(list.of.studies,".se"))]
                              require(reshape)
                              coefs = melt(coefs)
                              names(coefs) <- c("cell","study","coef")
                              ses = melt(ses)
                              data.long = cbind(coefs,ses[,"value"])
                              names(data.long)<-c("cell","study","coef","se")
                              res = split(data.long, f=data$ALSPAC.cell_type)
                              res = lapply(res, function(x) rma.uni(slab=x$study,yi=x$coef,sei=x$se,method="FE",weighted=TRUE))
                              res
                              }
extract.and.merge.meta.analysis <-function(meta.res,data){
                                  require(plyr)
                                  data.meta = ldply(lapply(meta.res, function(x)unlist(c(x$b[[1]],x[c("se","pval","QE","QEp","I2","H2")]))))
                                  colnames(data.meta)<-c("cell","coef.fe","se.fe","p.fe","q.fe","het.p.fe","i2.fe","h2.fe")
                                  data = merge(data,data.meta,by.x="ALSPAC.cell_type",by.y="cell",all.x=T)
                                  data
                                  }

meta.cell.results <- fixed.effects.meta.analysis(list.of.studies = studies, data = cell.results)
meta.cell.results.dataframe <- extract.and.merge.meta.analysis(meta.res = meta.cell.results, data = cell.results)

write.csv(meta.cell.results.dataframe,"cells.meta.results.childhood.csv")

# Draw forest plots
pdf("cells.forest.plots.childhood.pdf",width=10,height=6)
for(i in 1:length(meta.cell.results)){
par(mar=c(4,5,1,4))
forest(meta.cell.results[[i]],main="",digits=4,mlab="",xlab="Difference in cell proportion per 1SD increase in paternal BMI",xlim=c(-0.02,0.02),alim=c(-0.02,0.02),cex=1,ylim=c(-1,9))
text(-0.02,-1,pos=4,cex=1,
bquote(paste("FE Model (Q = ",
    .(formatC(meta.cell.results[[i]]$QE, digits=2, format="f")), ", df = ", .(meta.cell.results[[i]]$k - meta.cell.results[[i]]$p),
 ", p = ", .(formatC(meta.cell.results[[i]]$QEp, digits=2, format="f")), "; ", I^2, " = ",
  .(formatC(meta.cell.results[[i]]$I2, digits=1, format="f")), "%)")))
  text(0,8,paste0("\n\nMeta-analysis P-value = ",round(meta.cell.results[[i]]$pval,3)))
  text(0,8.2,cex=2,toupper(names(meta.cell.results))[i])

}
dev.off()

#Get Tau2 and I2 confidence intervals (only generated in RE meta-analysis):

Tau2.I2.confint.meta.analysis <- function(list.of.studies,data){
                              coefs = data[,c("ALSPAC.cell_type",paste0(list.of.studies,".effect"))]
							  colnames(coefs) <- c("cell",paste0(list.of.studies))
                              ses = data[,c("ALSPAC.cell_type",paste0(list.of.studies,".se"))]
                              require(reshape)
                              coefs = melt(coefs)
                              names(coefs) <- c("cell","study","coef")
                              ses = melt(ses)
                              data.long = cbind(coefs,ses[,"value"])
                              names(data.long)<-c("cell","study","coef","se")
                              res = split(data.long, f=data$ALSPAC.cell_type)
                              res = lapply(res, function(x) rma.uni(slab=x$study,yi=x$coef,sei=x$se,method="DL",weighted=TRUE))
                              res
                              }
Tau2.I2.confint <- lapply(Tau2.I2.confint.meta.analysis(list.of.studies = studies, data = cell.results),confint)
