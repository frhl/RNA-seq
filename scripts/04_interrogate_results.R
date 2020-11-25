# interrogating data results

# check results
files = list.files('derived/', pattern = '\\.csv', full.names = T)
lst = lapply(files, function(x) read.csv(x))
names(lst) <- basename(files)
lapply(lst, head)


# plot them!
library(ggplot2)
library(ggrepel)

df = read.csv("derived/201125_KO_SALIvsSALI.csv", row.names = NULL)
df$gene <- df$hgnc_symbol
df$label <- df$hgnc_symbol
df$label[!as.logical(abs(df$de_expression))] <- NA

## SALI vs KO sali
plt <- ggplot(df, aes(logFC, -log10(PValue), color = de_expression, label = label)) +
  ggtitle('SALI vs KO_SALI (FDR < 0.05)') +
  geom_point() +
  ggrepel::geom_label_repel(size = 2) +
  theme_bw()
plt




## check enrichment
library(genoppi)
pdf('201125_SALIvsKO_SALI_FDR005_GTEx_enrichment.pdf')
ref = gtex_table
df$significant <- df$FDR < 0.05
e <- calc_adjusted_enrichment(df, ref, intersectN = T)
colnames(e)[0] <- 'GTEx_tissue'
colnames(e)[9] <- 'FDR'
colnames(e)[10] <- 'FDR'




