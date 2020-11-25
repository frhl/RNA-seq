# run QC and check if
df <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
count.table <- df
rownames(count.table) <- df$ensembl_gene_id
count.table$hgnc_symbol <- NULL
count.table$ensembl_gene_id <- NULL

# cpm normalisation
normalised.count.table <- sweep(count.table, 2, colSums(count.table),FUN="/")
normalised.count.table <- normalised.count.table * 1000000

# Distance matrix heatmap
dists <- make.dist.heatmap(normalised.count.table)
plot.dist.heatmap(as.matrix(dists))

# PCA
PC.dat <- calc.pca(normalised.count.table)
plot.PCA.ids(PC.dat)



