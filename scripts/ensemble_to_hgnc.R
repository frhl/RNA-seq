#' @title esemble to HGNC mapping
#' @param

mapping <- read.csv('201124_fl_ensembl_to_hgnc_mapping.csv')[,2:3]
ensembl_to_hgnc <- function(x){
  res = unlist(lapply(x, function(x) mapping$hgnc_symbol[mapping$ensembl_gene_id == x]))
  res[nchar(res) < 1] <- NA
  return(res)
}
