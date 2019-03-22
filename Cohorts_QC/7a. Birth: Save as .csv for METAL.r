# Create dataframes for each model/cohort combination 
# (need two versions of every "mutual" dataframe, so that the number of models == number of dataframes in the list))

prepare.QCd.ewas.results <- function(result.file, ewas.dataframe, cohort.name){
key.result <- key[which(key$result==result.file),]
QCd.ewas.results <- ewas.dataframe[,intersect(colnames(ewas.dataframe), c(key.result$column.prefix,paste(key.result$column.prefix,key.result$column.suffix,sep=".")))]
colnames(QCd.ewas.results) <- c("n","coef","se","p","probeid")
filename <- paste0("meta/cohort_files_after_qc/birth.",cohort.name,".",result.file,".csv")
write.csv(QCd.ewas.results,filename,row.names=FALSE,quote=FALSE,eol="\r\n")
}

positions<-grep(names(All.EWAS$ALSPAC),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$ALSPAC),pattern="mutual"))))
ALSPAC2<- append(All.EWAS$ALSPAC,list(All.EWAS$ALSPAC[[positions[1]]]),positions[1])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[2]]]),positions[2])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[3]]]),positions[3])
ALSPAC2<- append(ALSPAC2,list(ALSPAC2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$BIB_asian),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$BIB_asian),pattern="mutual"))))
BIB_asian2<- append(All.EWAS$BIB_asian,list(All.EWAS$BIB_asian[[positions[1]]]),positions[1])
BIB_asian2<- append(BIB_asian2,list(BIB_asian2[[positions[2]]]),positions[2])
BIB_asian2<- append(BIB_asian2,list(BIB_asian2[[positions[3]]]),positions[3])

positions<-grep(names(All.EWAS$BIB_white),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$BIB_white),pattern="mutual"))))
BIB_white2<- append(All.EWAS$BIB_white,list(All.EWAS$BIB_white[[positions[1]]]),positions[1])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[2]]]),positions[2])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[3]]]),positions[3])
BIB_white2<- append(BIB_white2,list(BIB_white2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$CHAMACOS),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$CHAMACOS),pattern="mutual"))))
CHAMACOS2<- append(All.EWAS$CHAMACOS,list(All.EWAS$CHAMACOS[[positions[1]]]),positions[1])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[2]]]),positions[2])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[3]]]),positions[3])
CHAMACOS2<- append(CHAMACOS2,list(CHAMACOS2[[positions[4]]]),positions[4])

positions<-3
ENVIRONAGE2<- append(All.EWAS$ENVIRONAGE,list(All.EWAS$ENVIRONAGE[[positions[1]]]),positions[1])

positions<-grep(names(All.EWAS$GenerationR),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$GenerationR),pattern="mutual"))))
GenerationR2<- append(All.EWAS$GenerationR,list(All.EWAS$GenerationR[[positions[1]]]),positions[1])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[2]]]),positions[2])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[3]]]),positions[3])
GenerationR2<- append(GenerationR2,list(GenerationR2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$GOYA),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$GOYA),pattern="mutual"))))
GOYA2<- append(All.EWAS$GOYA,list(All.EWAS$GOYA[[positions[1]]]),positions[1])
GOYA2<- append(GOYA2,list(GOYA2[[positions[2]]]),positions[2])
GOYA2<- append(GOYA2,list(GOYA2[[positions[3]]]),positions[3])
GOYA2<- append(GOYA2,list(GOYA2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$INMA.nocombat),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$INMA.nocombat),pattern="mutual"))))
INMA.nocombat2<- append(All.EWAS$INMA.nocombat,list(All.EWAS$INMA.nocombat[[positions[1]]]),positions[1])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[2]]]),positions[2])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[3]]]),positions[3])
INMA.nocombat2<- append(INMA.nocombat2,list(INMA.nocombat2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$INMA.combat),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$INMA.combat),pattern="mutual"))))
INMA.combat2<- append(All.EWAS$INMA.combat,list(All.EWAS$INMA.combat[[positions[1]]]),positions[1])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[2]]]),positions[2])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[3]]]),positions[3])
INMA.combat2<- append(INMA.combat2,list(INMA.combat2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$MoBa1),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$MoBa1),pattern="mutual"))))
MoBa12<- append(All.EWAS$MoBa1,list(All.EWAS$MoBa1[[positions[1]]]),positions[1])
MoBa12<- append(MoBa12,list(MoBa12[[positions[2]]]),positions[2])
MoBa12<- append(MoBa12,list(MoBa12[[positions[3]]]),positions[3])
MoBa12<- append(MoBa12,list(MoBa12[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$MoBa2),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$MoBa2),pattern="mutual"))))
MoBa22<- append(All.EWAS$MoBa2,list(All.EWAS$MoBa2[[positions[1]]]),positions[1])
MoBa22<- append(MoBa22,list(MoBa22[[positions[2]]]),positions[2])
MoBa22<- append(MoBa22,list(MoBa22[[positions[3]]]),positions[3])
MoBa22<- append(MoBa22,list(MoBa22[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$MoBa3),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$MoBa3),pattern="mutual"))))
MoBa32<- append(All.EWAS$MoBa3,list(All.EWAS$MoBa3[[positions[1]]]),positions[1])
MoBa32<- append(MoBa32,list(MoBa32[[positions[2]]]),positions[2])
MoBa32<- append(MoBa32,list(MoBa32[[positions[3]]]),positions[3])
MoBa32<- append(MoBa32,list(MoBa32[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$PICCOLIPIU),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$PICCOLIPIU),pattern="mutual"))))
PICCOLIPIU2<- append(All.EWAS$PICCOLIPIU,list(All.EWAS$PICCOLIPIU[[positions[1]]]),positions[1])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[2]]]),positions[2])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[3]]]),positions[3])
PICCOLIPIU2<- append(PICCOLIPIU2,list(PICCOLIPIU2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$ProjectViva),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$ProjectViva),pattern="mutual"))))
ProjectViva2<- append(All.EWAS$ProjectViva,list(All.EWAS$ProjectViva[[positions[1]]]),positions[1])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[2]]]),positions[2])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[3]]]),positions[3])
ProjectViva2<- append(ProjectViva2,list(ProjectViva2[[positions[4]]]),positions[4])

positions<-grep(names(All.EWAS$RHEA),pattern="mutual")+c(0,seq(2:length(grep(names(All.EWAS$RHEA),pattern="mutual"))))
RHEA2<- append(All.EWAS$RHEA,list(All.EWAS$RHEA[[positions[1]]]),positions[1])
RHEA2<- append(RHEA2,list(RHEA2[[positions[2]]]),positions[2])
RHEA2<- append(RHEA2,list(RHEA2[[positions[3]]]),positions[3])
RHEA2<- append(RHEA2,list(RHEA2[[positions[4]]]),positions[4])

mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ALSPAC))]),ALSPAC2,"ALSPAC")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(BIB_asian))]),BIB_asian2,"BIB.asian")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(BIB_white))]),BIB_white2,"BIB.white")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(CHAMACOS))]),CHAMACOS2,"CHAMACOS")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ENVIRONAGE))]),ENVIRONAGE2,"ENVIRONAGE")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(GenerationR))]),GenerationR2,"GenerationR")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(GOYA))]),GOYA2,"GOYA")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(INMA.nocombat))]),INMA.nocombat2,"INMA.nocombat")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(INMA.combat))]),INMA.combat2,"INMA.combat")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(MoBa1))]),MoBa12,"MoBa1")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(MoBa2))]),MoBa22,"MoBa2")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(MoBa3))]),MoBa32,"MoBa3")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(PICCOLIPIU))]),PICCOLIPIU2,"PICCOLIPIU")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(ProjectViva))]),ProjectViva2,"ProjectViva")
mapply(prepare.QCd.ewas.results,unique(key$result[which(key$model %in% names(RHEA))]),RHEA2,"RHEA")
