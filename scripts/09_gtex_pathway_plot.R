
# load data
all.dirs <- c('derived/201126_H441_SUBvsKO_SUB',
          'derived/201126_H441_SUBvsSALI',
          'derived/201126_SALIvsKO_SALI',
          'derived/201126_KO_SUBvsKO_SALI')


pdf('derived/201126_GTEx_enrichment_only_lung.pdf', width = 8, height = 7)
for (p in all.dirs){
 
  # get name 
  dname = p
  raw_name = paste(unlist(lapply(strsplit(dname, '\\_'), function(x) x[2:length(x)])), collapse = '_')
  cond1 = unlist(lapply(strsplit(raw_name, 'vs'), function(x) x[1]))
  cond2 = unlist(lapply(strsplit(raw_name, 'vs'), function(x) x[2]))
  
  # compare all GTEx enrichment
  files = sort(list.files(p, pattern = 'GTEx', full.names = T))[1:2]
  pos = read.csv(files[grepl('pos', files)], sep = '\t')[,c(1,9)]
  pos$condition = as.factor(cond1)
  neg = read.csv(files[grepl('neg', files)], sep = '\t')[,c(1,9)]
  neg$condition = as.factor(cond2)
  
  #pos = pos[pos$list_name == 'Lung',]
  #neg = neg[neg$list_name == 'Lung',]
  
  # setup data.frame
  dt <- rbind(pos, neg)
  dt$tissue <- factor(dt$list_name, levels = unique(dt$list_name[order(dt$FDR)]))
  dt$list_name <- NULL
  dt$FDR <- -log10(dt$FDR)
  
  xmi <- -10
  xma <- 10
  
  p = ggplot(data = dt, aes(x = reorder(tissue, -FDR), fill = condition)) +
    geom_bar(stat = "identity", data = subset(dt, condition == cond1), aes(y=FDR)) +
    geom_bar(stat = "identity", data = subset(dt, condition == cond2), aes(y=FDR * (-1)) ) +
    scale_y_continuous(limits = c(xmi, xma), breaks = seq(xmi, xma, 10), labels = abs(seq(xmi, xma, 10))) + 
    theme(axis.text = element_text(colour = "black")) + 
    coord_flip() + 
    ggtitle(paste(raw_name,'- GTEx tissue enrichment'))+
    ylab("-log10(Hypergeometric FDR)") + 
    xlab("") +
    theme_minimal()
  print(p)
   
}
graphics.off()


