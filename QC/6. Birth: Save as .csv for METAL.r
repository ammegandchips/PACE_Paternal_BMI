# Create dataframes for each model/cohort combination 
# (need two versions of every "mutual" dataframe, so that the number of models == number of dataframes in the list))

prepare.QCd.ewas.results <- function(result.file, ewas.dataframe, cohort.name){
key.result <- key[which(key$result==result.file),]
QCd.ewas.results <- ewas.dataframe[,intersect(colnames(ewas.dataframe), c(key.result$column.prefix,paste(key.result$column.prefix,key.result$column.suffix,sep=".")))]
colnames(QCd.ewas.results) <- c("n","coef","se","p","probeid")
filename <- paste0("meta/cohort_files_after_qc/birth.",cohort.name,".",result.file,".csv")
write.csv(QCd.ewas.results,filename,row.names=FALSE,quote=FALSE,eol="\r\n")
}

positions<-grep(names(ALSPAC),pattern="mutual")+c(0,seq(2:length(grep(names(ALSPAC),pattern="mutual"))))
ALSPAC2<- append(ALSPAC,list(ALSPAC[[positions[1]]]),positions[1])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[2]]]),positions[2])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[3]]]),positions[3])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[4]]]),positions[4])

positions<-grep(names(BIB_asian),pattern="mutual")+c(0,seq(2:length(grep(names(BIB_asian),pattern="mutual"))))
BIB_asian2<- append(BIB_asian,list(BIB_asian[[positions[1]]]),positions[1])
BIB_asian2<- append(BIB_asian2,list(BIB_asian2[[positions[2]]]),positions[2])
BIB_asian2<- append(BIB_asian2,list(BIB_asian2[[positions[3]]]),positions[3])

positions<-grep(names(BIB_white),pattern="mutual")+c(0,seq(2:length(grep(names(BIB_white),pattern="mutual"))))
BIB_white2<- append(BIB_white,list(BIB_white[[positions[1]]]),positions[1])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[2]]]),positions[2])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[3]]]),positions[3])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[4]]]),positions[4])

positions<-grep(names(CHAMACOS),pattern="mutual")+c(0,seq(2:length(grep(names(CHAMACOS),pattern="mutual"))))
CHAMACOS2<- append(CHAMACOS,list(CHAMACOS[[positions[1]]]),positions[1])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[2]]]),positions[2])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[3]]]),positions[3])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[4]]]),positions[4])

positions<-3
ENVIRONAGE2<- append(ENVIRONAGE,list(ENVIRONAGE[[positions[1]]]),positions[1])

positions<-grep(names(GenerationR),pattern="mutual")+c(0,seq(2:length(grep(names(GenerationR),pattern="mutual"))))
GenerationR2<- append(GenerationR,list(GenerationR[[positions[1]]]),positions[1])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[2]]]),positions[2])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[3]]]),positions[3])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[4]]]),positions[4])

positions<-grep(names(INMA.nocombat),pattern="mutual")+c(0,seq(2:length(grep(names(INMA.nocombat),pattern="mutual"))))
INMA.nocombat2<- append(INMA.nocombat,list(INMA.nocombat[[positions[1]]]),positions[1])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[2]]]),positions[2])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[3]]]),positions[3])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[4]]]),positions[4])

positions<-grep(names(INMA.combat),pattern="mutual")+c(0,seq(2:length(grep(names(INMA.combat),pattern="mutual"))))
INMA.combat2<- append(INMA.combat,list(INMA.combat[[positions[1]]]),positions[1])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[2]]]),positions[2])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[3]]]),positions[3])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[4]]]),positions[4])

positions<-grep(names(PICCOLIPIU),pattern="mutual")+c(0,seq(2:length(grep(names(PICCOLIPIU),pattern="mutual"))))
PICCOLIPIU2<- append(PICCOLIPIU,list(PICCOLIPIU[[positions[1]]]),positions[1])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[2]]]),positions[2])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[3]]]),positions[3])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[4]]]),positions[4])

positions<-grep(names(ProjectViva),pattern="mutual")+c(0,seq(2:length(grep(names(ProjectViva),pattern="mutual"))))
ProjectViva2<- append(ProjectViva,list(ProjectViva[[positions[1]]]),positions[1])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[2]]]),positions[2])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[3]]]),positions[3])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[4]]]),positions[4])

positions<-grep(names(RHEA),pattern="mutual")+c(0,seq(2:length(grep(names(RHEA),pattern="mutual"))))
RHEA2<- append(RHEA,list(RHEA[[positions[1]]]),positions[1])
RHEA2<- append(RHEA2,list(RHEA2[[positions[2]]]),positions[2])
RHEA2<- append(RHEA2,list(RHEA2[[positions[3]]]),positions[3])
RHEA2<- append(RHEA2,list(RHEA2[[positions[4]]]),positions[4])

mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ALSPAC))]),ALSPAC2,"ALSPAC")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(BIB_asian))]),BIB_asian2,"BIB.asian")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(BIB_white))]),BIB_white2,"BIB.white")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(CHAMACOS))]),CHAMACOS2,"CHAMACOS")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ENVIRONAGE))]),ENVIRONAGE2,"ENVIRONAGE")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(GenerationR))]),GenerationR2,"GenerationR")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(INMA.nocombat))]),INMA.nocombat2,"INMA.nocombat")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(INMA.combat))]),INMA.combat2,"INMA.combat")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(PICCOLIPIU))]),PICCOLIPIU2,"PICCOLIPIU")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ProjectViva))]),ProjectViva2,"ProjectViva")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(RHEA))]),RHEA2,"RHEA")
