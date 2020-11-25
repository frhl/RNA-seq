# interrogating data results

# check results
files = list.files('derived/', pattern = '\\.csv', full.names = T)
lst = lapply(files, function(x) read.csv(x))
names(lst) <- basename(files)
lapply(lst, head)


# plot them!
library(ggplot2)
library(ggrepel)
library(genoppi)

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



pdf('SALIvsKO_SALI_')
ref = gtex_table
df$significant <- df$FDR < 0.05
e <- calc_adjusted_enrichment(df, ref, intersectN = T)
e$logFDR <- -log10(e$BH.FDR)
genoppi::plot_tissue_enrichment(e, 
                                col.tissue = 'list_name', 
                                col.value = 'logFDR',
                                xlab = 'GTEx (RNA tissue)',
                                ylab = '-log10(Hypergeometric FDR)')




