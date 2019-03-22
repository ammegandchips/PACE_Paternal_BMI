# Create dataframes for each model/cohort combination 
# (need two versions of every "mutual" dataframe, so that the number of models == number of dataframes in the list))

prepare.QCd.ewas.results <- function(result.file, ewas.dataframe, cohort.name){
key.result <- key[which(key$result==result.file),]
QCd.ewas.results <- ewas.dataframe[,intersect(colnames(ewas.dataframe), c(key.result$column.prefix,paste(key.result$column.prefix,key.result$column.suffix,sep=".")))]
colnames(QCd.ewas.results) <- c("n","coef","se","p","probeid")
filename <- paste0("meta/cohort_files_after_qc/childhood.",cohort.name,".",result.file,".csv")
write.csv(QCd.ewas.results,filename,row.names=FALSE,quote=FALSE,eol="\r\n")
}

positions<-grep(names(All.EWAS$ALSPAC),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$ALSPAC),pattern="mutual"))))
ALSPAC2<- append(All.EWAS$ALSPAC,list(All.EWAS$ALSPAC[[positions[1]]]),positions[1])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[2]]]),positions[2])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[3]]]),positions[3])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$CHAMACOS),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$CHAMACOS),pattern="mutual"))))
CHAMACOS2<- append(All.EWAS$CHAMACOS,list(All.EWAS$CHAMACOS[[positions[1]]]),positions[1])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[2]]]),positions[2])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[3]]]),positions[3])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$GenerationR),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$GenerationR),pattern="mutual"))))
GenerationR2<- append(All.EWAS$GenerationR,list(All.EWAS$GenerationR[[positions[1]]]),positions[1])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[2]]]),positions[2])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[3]]]),positions[3])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$HELIX),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$HELIX),pattern="mutual"))))
HELIX2<- append(All.EWAS$HELIX,list(All.EWAS$HELIX[[positions[1]]]),positions[1])
HELIX2<- append(HELIX2,list(HELIX2[[positions[2]]]),positions[2])
HELIX2<- append(HELIX2,list(HELIX2[[positions[3]]]),positions[3])
HELIX2<- append(HELIX2,list(HELIX2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$INMA),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$INMA),pattern="mutual"))))
INMA2<- append(All.EWAS$INMA,list(All.EWAS$INMA[[positions[1]]]),positions[1])
INMA2<- append(INMA2,list(INMA2[[positions[2]]]),positions[2])
INMA2<- append(INMA2,list(INMA2[[positions[3]]]),positions[3])
INMA2<- append(INMA2,list(INMA2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$ProjectViva),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$ProjectViva),pattern="mutual"))))
ProjectViva2<- append(All.EWAS$ProjectViva,list(All.EWAS$ProjectViva[[positions[1]]]),positions[1])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[2]]]),positions[2])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[3]]]),positions[3])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[4]]]),positions[4])

mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ALSPAC))]),ALSPAC2,"ALSPAC")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(CHAMACOS))]),CHAMACOS2,"CHAMACOS")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(GenerationR))]),GenerationR2,"GenerationR")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(INMA))]),INMA2,"INMA")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(HELIX))]),HELIX2,"HELIX")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ProjectViva))]),ProjectViva2,"ProjectViva")
