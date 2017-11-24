split.taxa = function (taxonomy) {
  require(stringr)

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
  
  return(splited)
  
}