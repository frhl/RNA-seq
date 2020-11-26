#setwd('project/')

# libs
library(edgeR)
library(ggplot2)

# ensemble to gene mapping
source('scripts/ensemble_to_hgnc.R')


# get matrix of raw expressions and ensemble IDs
df <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(df) <- df$ensembl_gene_id
df$hgnc_symbol <- NULL
df$ensembl_gene_id <- NULL

# conditions
layout <- read.table('data/layout_extended.txt', header = T)
conds <- factor(layout$Group)

# check data read depth
#read.depth <- apply(df, 2, sum)/1000000 
#read.depth <- colSums(df)/1000000 
#summary(read.depth)
#barplot(read.depth, las=2, cex.names=0.5, main = 'read depth (million transcrips)')

# check differentially expressed genes
y <- DGEList(counts=df, genes=row.names(df))

# library size?
quantile(y$samples$lib.size)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L,M) # quite small library size

# median libarary size divided by 10 
#keep.exprs <- filterByExpr(y, group=conds)
#y <- y[keep.exprs,, keep.lib.sizes=FALSE]
(cpm.threshold <- L/10)
keep <- rowSums(cpm(y)> cpm.threshold) >= 4
table(keep)
y <- y[keep,]

# normalized exprssion
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log=TRUE)
#boxplot(lcpm, las=2, col=col, main="")
#title(main="Normalised data",ylab="Log-cpm")

# setup design matrix
design <- model.matrix(~0+conds)
colnames(design) <- gsub("conds", "", colnames(design)) 
cont.matrix <- makeContrasts(H441_SUB-KO_SUB, # wt_sub vs ko_sub
                             H441_SUB-SALI,  # wt_sub vs wt_sali
                             SALI-KO_SALI,  # ko_sali vs wt_sal
                             levels=design)

pdf('201120_voom_plots.pdf')
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

# run over all samples
result <- lapply(1:ncol(cont.matrix), function(i){
  
  # get contrasts and calc logFC
  name = rownames(cont.matrix)
  contr_name = paste0(name[cont.matrix[,i] == 1], 'vs', name[cont.matrix[,i] == -1], collase = '')
  fit <- glmLRT(glmfit, contrast=cont.matrix[,i])
  de <- decideTestsDGE(fit) # FDR < 0.05
  
  # remove contaminants (keratins) in the samples
  fit.tabl <- fit$table

  # get FDR
  FDR <- p.adjust(fit.tabl$PValue, method="fdr")
  fit.tabl <- cbind(fit.tabl, FDR)
  fit.tabl$ensembl_gene_id <- rownames(fit.tabl)
  final <- merge(fit.tabl, mart, all.x = T)
  final <- final[,c(7,1,2:6)]
  
  # what genes are up/down regulated
  reg <- as.data.frame(de)
  reg$ensembl_gene_id <- rownames(reg)
  reg$de_expression <- reg[,1]
  reg[,1] <- NULL #reg$`1*H441_SUB -1*KO_SUB` <- NULL
  rownames(reg) <- NULL
  final <- merge(final, reg, all.x = T)
  
  # sort by FDR
  o <- order(final$FDR)
  final <- final[o,]
  
  # save to file
  outpath = paste0('derived/201125_',contr_name,'.csv')
  write(outpath, stdout())
  write.csv(final, outpath, row.names = F)
  
})
graphics.off()



