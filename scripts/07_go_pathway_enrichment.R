
library(ggplot2)
library(genoppi)


# get GO table
ref = goa_bp_table
colnames(ref) <- c('gene', 'ID', 'pathway')
ref = ref[,c(1,3)]
ref$significant <- TRUE
go_table <- ref

calc_go_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  e <- calc_adjusted_enrichment(df1, go_table, col.by = 'pathway', intersectN = F)
  colnames(e)[0] <- 'go_bp'
  colnames(e)[9] <- 'FDR'
  return(e)
}

# read files
files = list.files('derived/', pattern = '\\.csv', full.names = T)
lst = lapply(files, function(x) read.csv(x))
names(lst) <- basename(files)

for (f in files){
  
  # setup paths
  dname = tools::file_path_sans_ext(basename(f))
  dirpath = paste0('derived/',dname,'/')
  dir.create(dirpath)
  write(dname, stdout())
  
  # setup data
  df = read.csv(f, row.names = NULL)
  df$gene <- df$hgnc_symbol
  df$label <- df$hgnc_symbol
  df$label[!as.logical(abs(df$de_expression))] <- NA
  
  # calculat GO biological pathway enrichment
  df1 <- df
  df1$significant <- df1$FDR < 0.05
  both <- calc_go_enrichment(df1)
  df1$significant <- df1$FDR < 0.05 & df1$logFC > 0
  positive <- calc_go_enrichment(df1)
  df1$significant <- df1$FDR < 0.05 & df1$logFC < 0
  negative <- calc_go_enrichment(df1)
  
  ## check enrichment
  file = paste0(dirpath, "go_bp_enrichment_")
  write.table(both, paste0(file,"FDR_lt_005.tsv"), sep = '\t', quote = F)
  write.table(positive, paste0(file,"FDR_lt_005_logFC_pos.tsv"), sep = '\t', quote = F)
  write.table(negative, paste0(file,"FDR_lt_005_logFC_neg.tsv"), sep = '\t', quote = F)
  
  
}








