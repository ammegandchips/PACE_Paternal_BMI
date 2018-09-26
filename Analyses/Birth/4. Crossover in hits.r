#Are any hits for maternal/paternal BMI or males/females the same?

#pat vs mat
intersect(list.of.results$covs.pat[which(list.of.results$covs.pat$Pvalue<1e-5),"MarkerName"],list.of.results$covs.mat[which(list.of.results$covs.mat$Pvalue<1e-5),"MarkerName"])#0

#patmat vs matpat
intersect(list.of.results$covs.patmat[which(list.of.results$covs.patmat$Pvalue<1e-5),"MarkerName"],list.of.results$covs.matpat[which(list.of.results$covs.matpat$Pvalue<1e-5),"MarkerName"])#0

#males patmat vs matpat
intersect(list.of.results$boys.patmat[which(list.of.results$boys.patmat$Pvalue<1e-5),"MarkerName"],list.of.results$boys.matpat[which(list.of.results$boys.matpat$Pvalue<1e-5),"MarkerName"])#0

#females patmat vs matpat
intersect(list.of.results$girls.patmat[which(list.of.results$girls.patmat$Pvalue<1e-5),"MarkerName"],list.of.results$girls.matpat[which(list.of.results$girls.matpat$Pvalue<1e-5),"MarkerName"])#0

#patmat males vs females
patmat.girls.sites <- list.of.results$girls.patmat[which(list.of.results$girls.patmat$Pvalue<1e-5),]
patmat.boys.sites <- list.of.results$boys.patmat[which(list.of.results$boys.patmat$Pvalue<1e-5),]
intersect(patmat.girls.sites$MarkerName,patmat.boys.sites$MarkerName) 

#matpat males vs females
matpat.girls.sites <- list.of.results$girls.matpat[which(list.of.results$girls.matpat$Pvalue<1e-5),]
matpat.boys.sites <- list.of.results$boys.matpat[which(list.of.results$boys.matpat$Pvalue<1e-5),]
intersect(matpat.girls.sites$MarkerName,matpat.boys.sites$MarkerName) # 0
