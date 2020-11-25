
### helpers
source("http://www.well.ox.ac.uk/bioinformatics/training/RNASeq_Nov2020/scripts/simple_RNA_QC.R")

# run simple QC
layout.table <- read.table('data/layout_extended.txt')
count.table <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(count.table) <- df$ensembl_gene_id
count.table$hgnc_symbol <- NULL
count.table$ensembl_gene_id <- NULL

# GMS tool for QC analysis
simple.RNA.QC(count.table, layout.table, 'derived/201125_qc_summary')

# cpm normalisation
#normalised.count.table <- sweep(count.table, 2, colSums(count.table),FUN="/")
#normalised.count.table <- normalised.count.table * 1000000


