setwd('project/')

# biomart to generate mapping
library("biomaRt")
listMarts()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
attributes[grepl('hgnc', attributes$name),]
mart = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)

# apparantly, 10 and 16 have been switched around. Let's switch them back
layout <- read.csv('data/layout.txt', sep = '\t')
layout$Group[layout$Name == 10] <- "H441_SUB" # <--- changed only group
layout$Group[layout$Name == 16] <- "KO_SALI"  # <--- changed only group
layout$Name <- as.numeric(layout$Name)
layout$label <- paste0(layout$Group,'_rep',layout$Rep, sep = '')
layout <- layout[order(layout$Name), ]
length(unique(layout$label)) == 24 # expect 24 uniuqe column anames
write.table(layout, 'data/layout_extended.txt', sep = '\t')

# get gene count matrix / rename columns
df = read.csv('data/P200321_gene_count_table.txt', sep = '\t', header = T)
colnames(df)[1] <- "ensembl_gene_id"
colnames(df)[2:25] <- layout$label

# merge with biomart
dfmerge = merge(df, mart, all.x = T)
sum(is.na(dfmerge$hgnc_symbol))/nrow(dfmerge) # 7% not mapped
nrow(dfmerge)/nrow(df)
dfmerge <- dfmerge[,c(1,26,2:25)]
write.table(dfmerge, 'data/201124_gene_count_table_hgnc.txt', sep = '\t')




