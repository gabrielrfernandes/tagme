#' Split the taxonomic assinment into 7 taxa: Domain, Phylum, Class, Order, Family and Specie
#'
#' @param taxonomy The taxonomy object produced by tagmeFromFasta or tagmeFromDada2 functions.
#'
#' @return A data frame with 7 columns: Domain, Phylum, Class, Order, Family and Specie, describing the taxonmy for each sequence. This data frame can be used as input for phyloseq::tax_table().
#'
#' @export

split_taxa = function (taxonomy) {

  require("stringr")

  splited = stringr::str_split_fixed(taxonomy$Taxonomy, pattern = ";", 7)
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
