# interrogating data results

# libs
library(ggrepel)
library(genoppi) 
library(gplots)
library(ggplot2)
library(RColorBrewer)
source('scripts/enrichment_helpers.R')




date <- '201126'


# colors for plotting
cbp1 <- c("#E69F00", "#56B4E9", "#999999", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gms_colors <- c("X.NA" = "#999999", 
                "KO_SALI" = "#E69F00", 
                "SALI" = "#56B4E9",
                "H441_SUB" = "#009E73",
                "KO_SUB" = "#D55E00")

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
  
  # setup colors
  
  # label top 10 genes on each side
  dftop10 <- df
  dftop10$top <- F
  dftop10$top[df$de_expression == -1][1:15] <- T
  dftop10$top[df$de_expression == 1][1:15] <- T
  dftop10$label[dftop10$top == FALSE] <- NA
    
  
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
  
  # write conditions to files
  write.table(df_cond1, file = paste0(dirpath, date,'_',cond1,'_cond1','.tsv'), sep = '\t', quote = F, row.names = F)
  write.table(df_cond1_background, file = paste0(dirpath, date,'_',cond1,'_cond1_background','.tsv'), sep = '\t', quote = F, row.names = F)
  write.table(df_cond2, file = paste0(dirpath, date,'_',cond2,'_cond2','.tsv'), sep = '\t', quote = F, row.names = F)
  write.table(df_cond2_background, file = paste0(dirpath, date,'_',cond2,'_cond2_background','.tsv'), sep = '\t', quote = F, row.names = F)
  write.table(df_cond3, file = paste0(dirpath, date,'_cond3','.tsv'), sep = '\t', quote = F, row.names = F)
  write.table(df_cond3_background, file = paste0(dirpath, date,'_cond3_background','.tsv'), sep = '\t', quote = F, row.names = F)
  
  #if (F){
  
  # annotations (for plotting)
  annotations <- data.frame(
    xpos = c(-Inf, Inf),
    ypos =  c(Inf, Inf),
    annotateText = c(paste(n_cond2,'genes'), paste(n_cond1, 'genes')),
    color = c('black', 'blue'),
    hjustvar = c(-2,2) ,
    vjustvar = c(3,3))
  
  # open plotting device
  outplot <- paste0(dirpath,dname,'.pdf')
  pdf(outplot, width = 9, height = 7)
  
  # make plots
  plt1 <- ggplot(df, aes(logFC, -log10(PValue), color = condition, label = label)) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
    geom_point() +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    #scale_color_manual(values=cbp1) +
    ggrepel::geom_label_repel(size = 2, show.legend = FALSE) +
    scale_colour_manual(values = gms_colors) +
    theme_minimal()
  print(plt1)

  
  plt2 <- ggplot(dftop10, aes(logFC, -log10(PValue), color = condition, label = label)) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
    geom_point() +
    ggtitle(paste(dname,'(FDR < 0.01)')) +
    scale_colour_manual(values = gms_colors) +
    ggrepel::geom_label_repel(size = 2, show.legend = FALSE) +
    theme_minimal()
  print(plt2)
  
  #plt3 <- ggplot(df0001, aes(logFC, -log10(PValue), color = condition, label = label)) +
  #  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
  #  geom_point() +
  #  ggtitle(paste(dname,'(FDR < 0.01)')) +
  #  scale_color_manual(values=cbp1) +
  #  ggrepel::geom_label_repel(size = 2, show.legend = FALSE) +
  #  theme_minimal()
  #print(plt3)
  
  plt4 <- ggplot(df, aes(logFC, -log10(PValue), color = condition, label = label)) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
    geom_point() +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_colour_manual(values = gms_colors) +
    theme_minimal()
  print(plt4)
  
  ## only label things
  
  ## logFC x logCPM
  plt5 <- ggplot(df, aes(x = logCPM, y = logFC, color = condition, label = label)) +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_colour_manual(values = gms_colors) +
    geom_point() +
    theme_minimal()
  print(plt5)
  
  
  plt6 <- ggplot(dftop10, aes(x = logCPM, y = logFC, color = condition, label = label)) +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_colour_manual(values = gms_colors) +
    geom_point() +
    theme_minimal()
  print(plt6)
  
  # plot heatmaps of differentially expressed genes
  
  if (F){
  # calculate directional GTEx enrichment using a set
  # of hypergeometric tests with subsequent FDR correction.
  
  df1 <- df
  df1$significant <- df1$FDR < 0.05
  both <- calc_gtex_enrichment(df1)
  df1$significant <- df1$FDR < 0.05 & df1$logFC > 0
  positive <- calc_gtex_enrichment(df1)
  df1$significant <- df1$FDR < 0.05 & df1$logFC < 0
  negative <- calc_gtex_enrichment(df1)
  qq = calc_c5_enrichment(df1)
  

  
  ## check enrichment
  file = paste0(dirpath, "GTEx_enrichment_")
  write.table(both, paste0(file,"FDR_lt_005.tsv"), sep = '\t', quote = F)
  write.table(positive, paste0(file,"FDR_lt_005_logFC_pos.tsv"), sep = '\t', quote = F)
  write.table(negative, paste0(file,"FDR_lt_005_logFC_neg.tsv"), sep = '\t', quote = F)
  
  ## msigDB
  
  
  }
  
  graphics.off()
  
}








