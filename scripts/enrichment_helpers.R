### helpers start
calc_gtex_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  e <- calc_adjusted_enrichment(df, gtex_table, intersectN = T)
  colnames(e)[0] <- 'GTEx_tissue'
  colnames(e)[9] <- 'FDR'
  return(e)
}

### helpers start
calc_hpa_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  e <- calc_adjusted_enrichment(df, gtex_table, intersectN = T)
  colnames(e)[0] <- 'GTEx_tissue'
  colnames(e)[9] <- 'FDR'
  return(e)
}

### helpers start
calc_c5_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  ref = msigdb_c5_table
  colnames(ref) <- c('gene', 'pathway')
  ref$significant <- TRUE
  e <- calc_adjusted_enrichment(df, ref, col.by = 'gene', intersectN = F)
  colnames(e)[0] <- 'GTEx_tissue'
  colnames(e)[9] <- 'FDR'
  return(e)
}