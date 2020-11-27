
#
library(reshape2)
library(ggplot2)
library(ggpubr)



### helpers
source("http://www.well.ox.ac.uk/bioinformatics/training/RNASeq_Nov2020/scripts/simple_RNA_QC.R")

# run simple QC
layout.table <- read.table('data/layout_extended.txt')
count.table <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(count.table) <- df$ensembl_gene_id
count.table$hgnc_symbol <- NULL
count.table$ensembl_gene_id <- NULL

# GMS tool for QC analysis
simple.RNA.QC(count.table, layout.table, 'derived/201126_qc_summary')

# cpm normalisation
normalised.count.table <- sweep(count.table, 2, colSums(count.table),FUN="/")
normalised.count.table <- normalised.count.table * 1000000

heatmap(normalised.count.table)

# Look at surfactant gene expression
hgnc <- read.csv('201124_fl_ensembl_to_hgnc_mapping.csv')
hgnc$X <- NULL
count.table$ensembl_gene_id <- rownames(count.table)
colnames(hgnc); colnames(count.table)
count.table <- merge(count.table, hgnc, all.x = T)

# check for expression
mygenes = c('SFTPB', 'SFTPC', 'SFTPA1', 'SFTPA2', 'SFTPD')
d1 = count.table[count.table$hgnc_symbol %in% mygenes,]
rownames(d1) <- d1$hgnc_symbol
d1$ensembl_gene_id <- NULL
#d1$hgnc_symbol <- NULL
d1$ensemble_id <- NULL
d1 <- melt(d1, id.vars = 'hgnc_symbol')
colnames(d1) <- c('hgnc_symbol','condition', 'count')

# get new conditions names
nnames <- unlist(lapply(strsplit(as.character(d1$condition), '\\_'), function(x) paste(x[1:2], collapse = '_')))
nnames <- gsub('_rep[0-9]*','', nnames)
d1$cond <- nnames

# compare multiple surfactant genes
p <- ggplot(d1, aes(x=cond, y=count, fill = hgnc_symbol)) + 
  geom_boxplot() + coord_flip() + 
  xlab('Raw count') + ylab('Conditions') +
  ggtitle('Raw counts of Surfactant genes') +
  theme_minimal()
print(p)
ggsave(p, filename = 'derived/201125_raw_counts_of_surfactant_genes.pdf', width = 6, height = 10)

# compare target only
my_comparisons <- list( c("SALI", "KO_SUB"), c('SALI', 'KO_SALI'), c('SALI', 'H441_SUB'))
d2 <- d1[d1$hgnc_symbol == 'SFTPB',]
p2 <- ggplot(d2, aes(x=cond, y=count)) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  geom_point() +
  geom_boxplot() +
  xlab('Raw count') + ylab('Conditions') +
  ggtitle('Raw counts of SFTBP across cultures') +
  theme_minimal()
p2
ggsave(p2, filename = 'derived/201125_raw_counts_of_SFTBP.pdf', width = 6, height = 10)


# heatmap



