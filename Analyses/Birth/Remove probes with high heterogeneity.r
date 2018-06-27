Removing sites with high heterogeneity (HetPVal<1e-5)

identify.het <- function(ewas.dataframe){
  ewas.dataframe[which(ewas.dataframe$HetPVal<1e-3),"MarkerName"]
}

het.cpgs <- unique(unlist(lapply(list.of.results,identify.het)))#13424

filterprobes.het <- function(ewas.dataframe){
  ewas.dataframe <- ewas.dataframe[-which(ewas.dataframe$MarkerName %in% het.cpgs),]
  ewas.dataframe
}

list.of.results.het.removed <- lapply(list.of.results,filterprobes.het)

unlist(lapply(list.of.results,nrow))-unlist(lapply(list.of.results.het.removed,nrow))
