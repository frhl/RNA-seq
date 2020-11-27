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
  e <- calc_adjusted_enrichment(df, hpa_table, intersectN = T)
  colnames(e)[0] <- 'Human Protein Atlas'
  colnames(e)[9] <- 'FDR'
  return(e)
}

### helpers start
calc_h_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  ref = msigdb_h_table
  colnames(ref) <- c('gene', 'pathway')
  ref$significant <- TRUE
  e <- calc_adjusted_enrichment(df, ref, col.by = 'gene', intersectN = T)
  colnames(e)[0] <- ''
  colnames(e)[9] <- 'FDR'
  return(e)
}

### helpers start
calc_scRNA_enrichment <- function(df, scRNA){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  e <- calc_adjusted_enrichment(df, scRNA, col.by = 'cell', intersectN = T)
  colnames(e)[0] <- 'HPA_cell'
  colnames(e)[9] <- 'FDR'
  return(e)
}