# interrogating data results


library(genoppi) 
library(ggplot2)
library(RColorBrewer)


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
  
  ## compare contrast
  
  # count how many genes are expressed in either condition
  n_cond1 = sum(df$de_expression == 1)
  n_cond2 = sum(df$de_expression == -1)
  
  # annotations
  annotations <- data.frame(
    xpos = c(-Inf, Inf),
    ypos =  c(Inf, Inf),
    annotateText = c(paste(n_cond2,'genes'), paste(n_cond1, 'genes')),
    color = c('black', 'blue'),
    hjustvar = c(-2,2) ,
    vjustvar = c(3,3))
  
  # colors
  cbp1 <- c("#E69F00", "#56B4E9", "#999999", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # open plotting device
  outplot <- paste0(dirpath,dname,'.pdf')
  pdf(outplot, width = 12, height = 10)
  
  # make plots
  plt1 <- ggplot(df, aes(logFC, -log10(PValue), color = condition, label = label)) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
    geom_point() +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_color_manual(values=cbp1) +
    ggrepel::geom_label_repel(size = 1.5, show.legend = FALSE) +
    theme_minimal()
  print(plt1)
  
  plt2 <- ggplot(df, aes(logFC, -log10(PValue), color = condition, label = label)) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), color = 'black') +
    geom_point() +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_color_manual(values=cbp1) +
    theme_minimal()
  print(plt2)
  
  ## logFC x logCPM
  plt3 <- ggplot(df, aes(x = logCPM, y = logFC, color = condition, label = label)) +
    ggtitle(paste(dname,'(FDR < 0.05)')) +
    scale_color_manual(values=cbp1) +
    geom_point() +
    theme_minimal()
  print(plt3)
  
  # plot heatmaps of differentially expressed genes
  
  
  
  # calculate directional GTEx enrichment using a set
  # of hypergeometric tests with subsequent FDR correction.
  #df1 <- df
  #df1$significant <- df1$FDR < 0.05
  #both <- calc_gtex_enrichment(df1)
  #df1$significant <- df1$FDR < 0.05 & df1$logFC > 0
  #positive <- calc_gtex_enrichment(df1)
  #df1$significant <- df1$FDR < 0.05 & df1$logFC < 0
  #negative <- calc_gtex_enrichment(df1)
  
  ## check enrichment
  #file = paste0(dirpath, "GTEx_enrichment_")
  #write.table(both, paste0(file,"FDR_lt_005.tsv"), sep = '\t', quote = F)
  #write.table(positive, paste0(file,"FDR_lt_005_logFC_pos.tsv"), sep = '\t', quote = F)
  #write.table(negative, paste0(file,"FDR_lt_005_logFC_neg.tsv"), sep = '\t', quote = F)
  
  
  graphics.off()
  
}








