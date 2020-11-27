#setwd('project/')

# libs
library(edgeR)
library(ggplot2)

# ensemble to gene mapping
#source('scripts/ensemble_to_hgnc.R')

# for mapping genes
library("biomaRt")
listMarts()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
attributes[grepl('hgnc', attributes$name),]
mart = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)


# vars
date <- '201126'

################################
## setup design of experiment ##
################################

# this layout is just the sample 10 / 16 switched around.
layout <- read.table('data/layout_extended.txt', header = T)
conds <- factor(layout$Group)
design <- model.matrix(~0+conds)
colnames(design) <- gsub("conds", "", colnames(design)) 
cont.matrix <- makeContrasts(H441_SUB-KO_SUB, # wt_sub vs ko_sub
                             H441_SUB-SALI,  # wt_sub vs wt_sali
                             SALI-KO_SALI,  # ko_sali vs wt_sal
                             KO_SUB-KO_SALI,
                             levels=design)


#############################
## differential expression ##
#############################

# get matrix of raw expressions and ensemble IDs
df <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(df) <- df$ensembl_gene_id
df$hgnc_symbol <- NULL
df$ensembl_gene_id <- NULL

# Convert data into a DGEList with EdgeR
y <- DGEList(counts=df, genes=row.names(df))

# get an idea of the library size?
quantile(y$samples$lib.size)
size.mean <- mean(y$samples$lib.size) * 1e-6
size.median <- median(y$samples$lib.size) * 1e-6
print(size.mean, size.median)

# Set threshold to mean library size / 10
# in at least 4 samples (since we have 4 replicates)
(cpm.threshold <- size.mean/10)
keep <- rowSums(cpm(y)> cpm.threshold) >= 4
table(keep)
y <- y[keep,]
write.table(data.frame(keep), paste0('derived/',date, '_kept_rows.txt'))

# calculate norm factors
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log=TRUE)

# what are the top 10 % expressed genes
cpm(y, log = TRUE)

# Estimates a common negative binomial dispersion parameter for 
# a DGE dataset with a general experimental design.
par(mfrow=c(1,1))
y <- estimateGLMCommonDisp(y, design, verbose=T) # estimate of common dispersion
y <- estimateGLMTrendedDisp(y, design, verbose=T) # estimate of dispersion for each gene
y <- estimateGLMTagwiseDisp(y, design) # estimate of dispersion for each gene
plotBCV(y)

# summary of dispersions
disps <- list(y$common.dispersion, y$trended.dispersion, y$tagwise.dispersion)

# fit glm
glmfit <- glmFit(y, design)


#######################
## compare contrasts ##
#######################

# run over all samples
result <- lapply(1:ncol(cont.matrix), function(i){
  
  # get contrasts and calc logFC
  name = rownames(cont.matrix)
  contr_name = paste0(name[cont.matrix[,i] == 1], 'vs', name[cont.matrix[,i] == -1], collase = '')
  fit <- glmLRT(glmfit, contrast=cont.matrix[,i])
  de <- decideTestsDGE(fit) # FDR < 0.05
  table(de)
  
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
  outpath = paste0('derived/',date,'_',contr_name,'.csv')
  write(outpath, stdout())
  write.csv(final, outpath, row.names = F)
  
})
graphics.off()



