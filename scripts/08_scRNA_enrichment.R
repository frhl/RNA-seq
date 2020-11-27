





# interrogating data results

# libs
library(ggrepel)
library(genoppi) 
library(gplots)
library(ggplot2)
library(RColorBrewer)
source('scripts/enrichment_helpers.R')


###################
# prep scRNA data #
###################

if (F){
  
  scRNA <- read.csv('data/extdata/hpa_rna_single_cell_type.txt', sep = '\t')
  scRNA.cells <- unique(scRNA$Cell.type)
  scRNA.sig <-lapply(scRNA.cells, function(cell){
    print(cell)
    subdf <- scRNA[scRNA$Cell.type == cell,]
    subdf$significant <- subdf$NX >= quantile(subdf$NX, 0.9)
    return(subdf)
  })
  scRNA <- do.call(rbind, scRNA.sig)
  colnames(scRNA) <- c('ensemble_gene_id', 'gene', 'cell', 'NX', 'significant')
  write.table(scRNA, 'data/extdata/hpa_rna_single_cell_type_sig.txt', sep = '\t', quote = F)
  
}




date <- '201127'


# load results from file
files = list.files('derived/', pattern = '201126.+\\.csv', full.names = T, include.dirs = F, recursive = F)
lst = lapply(files, function(x) read.csv(x))
names(lst) <- basename(files)
lapply(lst, head)

#################################
# for each file, get some stats #
# ###############################

for (f in files){
  
  # setup paths
  dname = tools::file_path_sans_ext(basename(f))
  dirpath = paste0('derived/',dname,'/')
  dir.create(dirpath)
  write(dname, stdout())
  
  # setup named conditions (for plotting what side the differential gene expression belongs to)
  raw_name = paste(unlist(lapply(strsplit(dname, '\\_'), function(x) x[2:length(x)])), collapse = '_')
  cond1 = unlist(lapply(strsplit(raw_name, 'vs'), function(x) x[1]))
  cond2 = unlist(lapply(strsplit(raw_name, 'vs'), function(x) x[2]))
  
  # setup data
  df = read.csv(f, row.names = NULL)
  df$gene <- df$hgnc_symbol
  df$label <- df$hgnc_symbol
  df$label[!as.logical(abs(df$de_expression))] <- NA
  df$condition = 'X.NA'
  df$condition[df$de_expression == 1] <- cond1
  df$condition[df$de_expression == -1] <- cond2
  df$condition <- as.factor(df$condition)
  
  ######################
  ## compare contrast ##
  ######################
  
  # count how many genes are expressed in either condition
  n_cond1 = sum(df$de_expression == 1)
  n_cond2 = sum(df$de_expression == -1)
  
  # get df representing each condition for each condition
  df_cond1 = df[df$de_expression == 1,]$ensembl_gene_id
  df_cond1_background = df[df$de_expression != 1,]$ensembl_gene_id
  df_cond2 = df[df$de_expression == -1,]$ensembl_gene_id
  df_cond2_background = df[df$de_expression != -1,]$ensembl_gene_id 
  df_cond3 = df[abs(df$de_expression) == 1,]$ensembl_gene_id # FDR < 0.05, abs(logFC) > 0
  df_cond3_background = df[abs(df$de_expression) != 1,]$ensembl_gene_id # FDR < 0.05, abs(logFC) > 0
  
  if (T){
    
    df1 <- df
    df1$significant <- df1$FDR < 0.05
    both <- calc_scRNA_enrichment(df1, scRNA)
    
    df1$significant <- df1$FDR < 0.05 & df1$logFC > 0
    positive <- calc_scRNA_enrichment(df1, scRNA)
    
    df1$significant <- df1$FDR < 0.05 & df1$logFC < 0
    negative <- calc_scRNA_enrichment(df1, scRNA)
    
    ## scRNA HPA enrichment
    scfile = paste0(dirpath, "scRNA_HPA_enrichment_")
    write.table(both, paste0(scfile,"FDR_lt_005.tsv"), sep = '\t', quote = F)
    write.table(positive, paste0(scfile,"FDR_lt_005_logFC_pos.tsv"), sep = '\t', quote = F)
    write.table(negative, paste0(scfile,"FDR_lt_005_logFC_neg.tsv"), sep = '\t', quote = F)
    
  }
  
  graphics.off()
  
}








