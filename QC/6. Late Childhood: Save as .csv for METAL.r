# Create dataframes for each model/cohort combination 
# (need two versions of every "mutual" dataframe, so that the number of models == number of dataframes in the list))

prepare.QCd.ewas.results <- function(result.file, ewas.dataframe, cohort.name){
key.result <- key[which(key$result==result.file),]
QCd.ewas.results <- ewas.dataframe[,intersect(colnames(ewas.dataframe), c(key.result$column.prefix,paste(key.result$column.prefix,key.result$column.suffix,sep=".")))]
colnames(QCd.ewas.results) <- c("n","coef","se","p","probeid")
filename <- paste0("meta/cohort_files_after_qc/late_childhood.",cohort.name,".",result.file,".csv")
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

mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(All.EWAS$ALSPAC))]),ALSPAC2,"ALSPAC")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(All.EWAS$CHAMACOS))]),CHAMACOS2,"CHAMACOS")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(All.EWAS$GenerationR))]),GenerationR2,"GenerationR")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(All.EWAS$HELIX))]),HELIX2,"HELIX")
