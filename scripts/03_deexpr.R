setwd('project/')

# libs
library(edgeR)
library(ggplot2)

# ensemble to gene mapping
mapping <- read.csv('201124_fl_ensembl_to_hgnc_mapping.csv')[,2:3]
ensembl_to_hgnc <- function(x){
  res = unlist(lapply(x, function(x) mapping$hgnc_symbol[mapping$ensembl_gene_id == x]))
  res[nchar(res) < 1] <- NA
  return(res)
}


### pre-processing

# get matrix of raw expressions and ensemble IDs
df <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(df) <- df$ensembl_gene_id
df$hgnc_symbol <- NULL
df$ensembl_gene_id <- NULL

# conditions
layout <- read.table('data/layout_extended.txt', header = T)
conds <- factor(layout$Group)

# check data read depth
read.depth <- apply(df, 2, sum)/1000000 
read.depth <- colSums(df)/1000000 
summary(read.depth)
barplot(read.depth, las=2, cex.names=0.5, main = 'read depth (million transcrips)')

# check differentially expressed genes
y <- DGEList(counts=df, genes=row.names(df))

# library size?
quantile(y$samples$lib.size)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L,M) # quite small library size

# keep gene rows if there at least 2 experiments (~25%) with at least 200.000 counts
keep.exprs <- filterByExpr(y, group=conds)
#y <- y[keep.exprs,, keep.lib.sizes=FALSE]
(cpm.threshold <- L/10)
keep <- rowSums(cpm(y)> cpm.threshold) >= 10
table(keep)
y <- y[keep,]

# unnormalized expression
par(mfrow=c(1,2))
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")

# normalized exprssion
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

# setup design matrix
design <- model.matrix(~0+conds)
colnames(design) <- gsub("conds", "", colnames(design)) 
cont.matrix <- makeContrasts(H441_SUB-KO_SUB, 
                             H441_SUB-SALI,
                             KO_SALI-SALI,
                             levels=design)
cont.matrix

# visual check on the level of filtering performed upstream
par(mfrow=c(1,2))
v <- voom(y, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
# no drop in variance observed in the lower end of the 
# average expression. Seems to be a fair threshold.

# what genes are differentially expressed under each condition?
summary(decideTests(efit))

# Estimates a common negative binomial dispersion parameter for 
# a DGE dataset with a general experimental design.
par(mfrow=c(1,1))
y <- estimateGLMCommonDisp(y, design, verbose=T) 
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

# fit glm
glmfit <- glmFit(y, design)

# look at first differential expression
H441_SUBvsKO_sub <- glmLRT(glmfit, contrast=cont.matrix[,1])
de <- decideTestsDGE(H441_SUBvsKO_sub)

# what genes are upregulated / downregulated?
summary(de <- decideTestsDGE(H441_SUBvsKO_sub, adjust.method = "fdr"))
detags <- rownames(y)[as.logical(de)]
(detags_hgnc <- sort(ensembl_to_hgnc(detags)))
'SFTPB' %in% detags_hgnc # present!

plotSmear(H441_SUBvsKO_sub, de.tags=detags, main="H441_SUBvsKO_sub", ylim=c(-10,10))
abline(h=c(-1,1), col="blue")

# order table by pvalue
o <- order(H441_SUBvsKO_sub$table$PValue)
H441_SUBvsKO_sub <- H441_SUBvsKO_sub$table[o,]

## compute and add adjusted p-values
adjp <- p.adjust(H441_SUBvsKO_sub$PValue, method="fdr")
H441_SUBvsKO_sub <- cbind(H441_SUBvsKO_sub, adjp)
H441_SUBvsKO_sub$ensembl_gene_id <- rownames(H441_SUBvsKO_sub)

# merge mapping
H441_SUBvsKO_sub <- merge(H441_SUBvsKO_sub, mapping)
ggplot(H441_SUBvsKO_sub, aes(x = logFC, y = -log10(PValue))) +
  geom_point()


### look at different coditions

result <- lapply(1:ncol(cont.matrix), function(i){
  
  contr <- paste(rownames(cont.matrix)[as.logical(abs(cont.matrix[,i]))], collapse = 'vs')
  fit <- glmLRT(glmfit, contrast=cont.matrix[,i])
  de <- decideTestsDGE(fit)
  
  o <- order(fit$table$PValue)
  fit.tabl <- fit$table[o,]
  
  # 
  FDR <- p.adjust(fit.tabl$PValue, method="fdr")
  fit.tabl <- cbind(fit.tabl, FDR)
  fit.tabl$ensembl_gene_id <- rownames(fit.tabl)
  fit.tabl$hgnc <- ensembl_to_hgnc(fit.tabl$ensembl_gene_id)
  fit.tabl[]
  
  #mart = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters = values = fit.tabl$ensembl_gene_id, mart = ensembl)
  
})




