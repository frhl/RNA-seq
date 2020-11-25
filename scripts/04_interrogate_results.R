# interrogating data results

library(genoppi) # for enrichment analysis
library(ggplot2)


### helpers start

calc_gtex_enrichment <- function(df){
  # hypergeometric overlap analysis adjusted with FDR
  # data from genoppi package library(genoppi)
  e <- calc_adjusted_enrichment(df, gtex_table, intersectN = T)
  colnames(e)[0] <- 'GTEx_tissue'
  colnames(e)[9] <- 'FDR'
  return(e)
}

### helpers end

# check results
files = list.files('derived/', pattern = '\\.csv', full.names = T)
lst = lapply(files, function(x) read.csv(x))
names(lst) <- basename(files)
lapply(lst, head)

# plot them!
library(ggplot2)
library(ggrepel)


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
  
  ## compare contrast
  
  ## logFC x -log10(pvalue)
  outplot <- paste0(dirpath,dname,'.pdf')
  pdf(outplot, width = 12, height = 14)
  plt <- ggplot(df, aes(logFC, -log10(PValue), color = de_expression, label = label)) +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    geom_point() +
    ggrepel::geom_label_repel(size = 2) +
    theme_bw()
  print(plt)
  
  ## logFC x logCPM
  plt <- ggplot(df, aes(x = logCPM, y = logFC, color = de_expression, label = label)) +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    geom_point() +
    ggrepel::geom_label_repel(size = 2) +
    theme_bw()
  print(plt)
  
  # calc directional GTEx enrichment using a set
  # of hypergeometric tests with subsequent FDR correction.
  df1 <- df
  df1$significant <- df1$FDR < 0.05
  both <- calc_gtex_enrichment(df1)
  #genoppi::plot_tissue_enrichment(data, col.tissue = 'listName', col.value = '')
  #df1$significant <- df1$FDR < 0.05 & df1$logFC > 0
  #positive <- calc_gtex_enrichment(df1)
  #df1$significant <- df1$FDR < 0.05 & df1$logFC < 0
  #negative <- calc_gtex_enrichment(df1)
  
  ## check enrichment
  file = paste0(dirpath, "GTEx_enrichment_")
  write.csv(both, paste0(file,"FDR_lt_005.tsv"), sep = '\t', quote = F)
  #write.csv(positive, paste0(file,"FDR_lt_005_pos.tsv"), sep = '\t', quote = F)
  #write.csv(negative, paste0(file,"FDR_lt_005_neg.tsv"), sep = '\t', quote = F)
  
  graphics.off()
  
}








