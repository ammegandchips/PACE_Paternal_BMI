IGF2 <- annotation[which(annotation$chromosome=="chr11" & (annotation$position>2150342 & annotation$position<2170833)),]
MEST <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>130126012 & annotation$position<130146133)),]
PEG3 <- annotation[which(annotation$chromosome=="chr19" & (annotation$position>57321445 & annotation$position<57352096)),]
NNAT <- annotation[which(annotation$chromosome=="chr20" & (annotation$position>36149607 & annotation$position<36152092)),]
NDN <- annotation[which(annotation$chromosome=="chr15" & (annotation$position>23930554 & annotation$position<23932450)),]
SNRPN <- annotation[which(annotation$chromosome=="chr15" & (annotation$position>25068794 & annotation$position<25664609)),]
SGCE <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>94214536 & annotation$position<94285521)),]
PEG10 <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>94285637 & annotation$position<94299007)),]
MEG3 <- annotation[which(annotation$chromosome=="chr14" & (annotation$position>101245747 & annotation$position<101327368)),]
H19 <- annotation[which(annotation$chromosome=="chr11" & (annotation$position>2016406 & annotation$position<2022700)),]
PLAGL1 <- annotation[which(annotation$chromosome=="chr6" & (annotation$position>144261437 & annotation$position<144385735)),]
GRB10 <- annotation[which(annotation$chromosome=="chr7" & (annotation$position>50657760 & annotation$position<50861159)),]

list.of.results$covs.patmat <- list.of.results$covs.patmat[order(list.of.results$covs.patmat$Effect),]

imprinted.regions <-  rbind(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),],
list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),])

imprinted.regions$imprinted.region <- c(rep("IGF2",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),])),
rep("MEST",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),])),
rep("PEG3",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),])),
rep("NNAT",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),])),
rep("NDN",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),])),
rep("SNRPN",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),])),
rep("SGCE",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),])),
rep("PEG10",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),])),
rep("MEG3",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),])),
rep("H19",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),])),
rep("PLAGL1",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),])),
rep("GRB10",nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),]))
)

imprinted.regions$ID <- c(1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% IGF2$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEST$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG3$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NNAT$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% NDN$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SNRPN$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% SGCE$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PEG10$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% MEG3$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% H19$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% PLAGL1$name),]),
1:nrow(list.of.results$covs.patmat[which(list.of.results$covs.patmat$MarkerName %in% GRB10$name),]))


require(ggplot2)

Plot<-ggplot(imprinted.regions, aes(x=1, y=ID))+
geom_tile(aes(fill = Effect))+
facet_wrap(~imprinted.region, ncol=3,scales="free")+
xlab("")+ylab("")+
scale_fill_gradientn(colours = c("darkblue", "white", "red"),name="Effect\nestimate",limits=c(-max(abs(imprinted.regions$Effect)),max(abs(imprinted.regions$Effect))))+
theme_grey(base_size=8) + labs(x = "",y = "") + 
scale_x_discrete(expand = c(0, 0)) +
scale_y_discrete(expand = c(0, 0)) + theme(axis.ticks = element_blank(),axis.text.y = element_text(colour="black",size=11))+
theme(strip.background = element_rect(fill="white"), strip.text=element_text(size=12))+
ggtitle("Associations between paternal BMI and \noffspring cord blood methylation at imprinted regions\n")+
theme(plot.title=element_text(size=16))+
theme(legend.text = element_text(size = 12),legend.title=element_text(size=12))
Plot

pdf("imprinted.patmatbmi.pdf")
Plot
dev.off()
