# sensitivity analysis: how many genes are found to be upregulated
# with respect to different thresholds of CPM and columns?

# get matrix of raw expressions and ensemble IDs
df <- read.table('data/201124_gene_count_table_hgnc.txt', sep = '\t')
rownames(df) <- df$ensembl_gene_id
df$hgnc_symbol <- NULL
df$ensembl_gene_id <- NULL

# conditions
layout <- read.table('data/layout_extended.txt', header = T)
conds <- factor(layout$Group)

# thresholds
thresholds = seq(0.2, 1, by = 0.05)
cols = 1:15
lst <- list()


for (t in thresholds){
  
  for (q in cols){
    
    # check differentially expressed genes
    y <- DGEList(counts=df, genes=row.names(df))
    
    # library size?
    quantile(y$samples$lib.size)
    L <- mean(y$samples$lib.size) * 1e-6
    M <- median(y$samples$lib.size) * 1e-6
    c(L,M) # quite small library size
    
    #y <- y[keep.exprs,, keep.lib.sizes=FALSE]
    #(cpm.threshold <- L/10)
    cpm.threshold = t
    keep <- rowSums(cpm(y)> cpm.threshold) >= q
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
    
    #pdf('201120_voom_plots.pdf')
    # visual check on the level of filtering performed upstream
    par(mfrow=c(1,2))
    v <- voom(y, design, plot=TRUE)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=cont.matrix)
    efit <- eBayes(vfit)
    
    # get doown/upregulated
    mat <- summary(decideTests(efit))
    print(mat)
    
    m1 <- data.frame(cpm = t, rows = q, down = mat[1,1], NotSig = mat[2,1], up = mat[3,1])
    m2 <- data.frame(cpm = t, rows = q, down = mat[1,2], NotSig = mat[2,2], up = mat[3,2])
    m3 <- data.frame(cpm = t, rows = q, down = mat[1,3], NotSig = mat[2,3], up = mat[3,3])
    
    lst[['H441_SUB - KO_SUB']] <- rbind(lst[['H441_SUB - KO_SUB']], m1)
    lst[['H441_SUB - SALI']] <- rbind(lst[['H441_SUB - SALI']], m2)
    lst[['SALI - KO_SALI']] <- rbind(lst[['SALI - KO_SALI']], m3)

  }
  
}

# create dirs
dir.create('derived/sensitvity/')
lapply(names(lst), function(x) {
  write.csv(lst[[x]], paste0('derived/sensitvity/201125_sensitivity_',gsub(' ','',x), '.csv'))
})



#mat <- do.call(rbind, lst)

graphics.off()
pdf('derived/sensitvity/201125_sensitivity_plots.pdf', width =  12, height = 8)
for (i in 1:3){
  
  d1 <- lst[[i]]
  d1$NotSig <- NULL
  d1 <- melt(d1, id.vars = c('cpm', 'rows'))
  plt1 <- ggplot(d1, aes(x = cpm, y = rows, fill = value, label = value)) +
    geom_tile(color = 'black') +
    ggtitle(names(lst)[[i]]) +
    geom_text() +
    theme_minimal() +
    facet_grid(~variable)
  
  print(plt1)
}
graphics.off()





