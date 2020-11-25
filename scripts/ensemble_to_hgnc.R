#' @title esemble to HGNC mapping
#' @param



## FYI this is super slow, but works
mapping <- read.csv('201124_fl_ensembl_to_hgnc_mapping.csv')[,2:3]
ensembl_to_hgnc <- function(x){
  res = unlist(lapply(x, function(x) mapping$hgnc_symbol[mapping$ensembl_gene_id == x]))
  res[nchar(res) < 1] <- NA
  return(res)
}

## this is fast
mapping <- read.csv('201124_fl_ensembl_to_hgnc_mapping.csv')[,2:3]
df = read.csv('data/P200321_gene_count_table.txt', sep = '\t', header = T)
colnames(df)[1] <- "ensembl_gene_id"
colnames(df)[2:25] <- layout$label

# merge with biomart
dfmerge = merge(df, mapping, all.x = T)





