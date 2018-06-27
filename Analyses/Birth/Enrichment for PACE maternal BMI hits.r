#Enrichment for maternal BMI sites (74 of 86 available in this analysis)

myQQ <- function(p.sva, ci = 0.95,Title) {
  n  <- length(p.sva)
  df <- data.frame(
    observed = -log10(p.sva),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
	ylim(0,6) +
    geom_line(aes(expected, clower), linetype = 2) +
    ggtitle(Title) +
	theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(log10Pe) +
    ylab(log10Po)  
} 

previous86 <-read.csv("/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/CpGs_86.csv",stringsAsFactors=FALSE)
summary(previous86 %in% list.of.results.het.removed$covs.pat$MarkerName)
png("EnrichmentQQ.86.patbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.pat$Pvalue[which(list.of.results.het.removed$covs.pat$MarkerName %in% previous86$CpG)]),Title="Paternal BMI P-values:\nEnrichment for Sharp 2017 maternal BMI CpGs")
dev.off()
png("EnrichmentQQ.86.matbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.mat$Pvalue[which(list.of.results.het.removed$covs.mat$MarkerName %in% previous86$CpG)]),Title="Maternal BMI P-values:\nEnrichment for Sharp 2017 maternal BMI CpGs")
dev.off()

png("EnrichmentQQ.86.patmatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.patmat$Pvalue[which(list.of.results.het.removed$covs.patmat$MarkerName %in% previous86$CpG)]),Title="Paternal BMI (adjusted for maternal BMI) P-values:\nEnrichment for Sharp 2017 maternal BMI CpGs")
dev.off()
png("EnrichmentQQ.86.matpatbmi.png",width=500,height=500)
myQQ(p.sva = sort(list.of.results.het.removed$covs.matpat$Pvalue[which(list.of.results.het.removed$covs.matpat$MarkerName %in% previous86$CpG)]),Title="Maternal BMI (adjusted for paternal BMI) P-values:\nEnrichment for Sharp 2017 maternal BMI CpGs")
dev.off()
