#Are any hits for maternal/paternal BMI the same?

#pat vs mat
intersect(list.of.results.het.removed$covs.pat[which(list.of.results.het.removed$covs.pat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$covs.mat[which(list.of.results.het.removed$covs.mat$Pvalue<1e-5),"MarkerName"])#0

#patmat vs matpat
intersect(list.of.results.het.removed$covs.patmat[which(list.of.results.het.removed$covs.patmat$Pvalue<1e-5),"MarkerName"],list.of.results.het.removed$covs.matpat[which(list.of.results.het.removed$covs.matpat$Pvalue<1e-5),"MarkerName"])#0


