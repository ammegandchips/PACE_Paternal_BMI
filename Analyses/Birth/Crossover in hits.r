#Are any hits for maternal/paternal BMI the same?

#pat vs mat
intersect(list.of.results.het.removed$covs.pat[which(list.of.results.het.removed$covs.pat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$covs.mat[which(list.of.results.het.removed$covs.mat$Pvalue<1e-5),"MarkerName"])#0

#patmat vs matpat
intersect(list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$Pvalue<1e-5),"MarkerName"])#0

#males patmat vs matpat
intersect(list.of.results.het.removed$boys.patmat[which(list.of.results.het.removed$boys.patmat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$boys.matpat[which(list.of.results.het.removed$boys.matpat$Pvalue<1e-5),"MarkerName"])#0

#females patmat vs matpat
intersect(list.of.results.het.removed$girls.patmat[which(list.of.results.het.removed$girls.patmat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$girls.matpat[which(list.of.results.het.removed$girls.matpat$Pvalue<1e-5),"MarkerName"])#0

#patmat males vs females
patmat.girls.sites <- list.of.results.het.removed$girls.patmat[which(list.of.results.het.removed$girls.patmat$Pvalue<1e-5),]
patmat.boys.sites <- list.of.results.het.removed$boys.patmat[which(list.of.results.het.removed$boys.patmat$Pvalue<1e-5),]
intersect(patmat.girls.sites$MarkerName,patmat.boys.sites$MarkerName) 

#matpat males vs females
matpat.girls.sites <- list.of.results.het.removed$girls.matpat[which(list.of.results.het.removed$girls.matpat$Pvalue<1e-5),]
matpat.boys.sites <- list.of.results.het.removed$boys.matpat[which(list.of.results.het.removed$boys.matpat$Pvalue<1e-5),]
intersect(matpat.girls.sites$MarkerName,matpat.boys.sites$MarkerName) # 0
