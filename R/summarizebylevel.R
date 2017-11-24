summarize.by.level = function (table, level = "Genus", taxa.are.rows = T) {
  require(stringr)
  
  if (isTRUE(taxa.are.rows)) {
    rownames(table) = paste0("sq", seq(length(rownames(table))))
    table = data.frame(table)
  } else {
    colnames(table) = paste0("sq", seq(length(colnames(table))))
    table = data.frame(t(table))
  }
  
  splited = str_split_fixed(taxonomy$Taxonomy, pattern = ";", 7)
  splited = gsub(".__", "", splited)
  for (i in 1:dim(splited)[1]) {
    for (j in 1:7) {
      if (splited[i,j] != "") {
        last_assigned = splited[i,j]
      }
      else {
        splited[i,j] = paste("unclassified", last_assigned, sep = "_")
      }
    }
  }
  colnames(splited) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie")
  rownames(splited) = taxonomy$seqID
  
  table_aggregate = aggregate(table, by=list(as.factor(splited[,level])), FUN=sum)
  row.names(table_aggregate) = table_aggregate[,1]
  table_aggregate = table_aggregate[,-1]
  
  return(table_aggregate)
  
}