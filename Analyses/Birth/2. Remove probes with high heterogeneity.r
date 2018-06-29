#Removing sites with high heterogeneity (HetPVal<1e-5)

filterprobes.het <- function(ewas.dataframe){
  ewas.dataframe[which(ewas.dataframe$HetPVal>1e-3),]
  ewas.dataframe
}

list.of.results.het.removed <- lapply(list.of.results,filterprobes.het)

unlist(lapply(list.of.results,nrow))-unlist(lapply(list.of.results.het.removed,nrow))
