
expr <- read.table(
  "data_raw/GSE212406_gene_fpkm_matrix.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
file.info("data_raw/GSE212406_gene_fpkm_matrix.txt")
dim(expr)
colnames(expr)
asthma <- expr[, c(
  "gene_id",
  "Control72h1",
  "Control72h2",
  "Control72h3",
  "OVA72h1",
  "OVA72h2",
  "OVA72h3",
  "OVA72h4"
)]
dim(asthma)
write.csv(
  asthma,
  "data_processed/Control72h_vs_OVA72h.csv",
  row.names = FALSE
)
metadata <- data.frame(
  Sample = c(
    "Control72h1",
    "Control72h2",
    "Control72h3",
    "OVA72h1",
    "OVA72h2",
    "OVA72h3",
    "OVA72h4"
  ),
  Condition = c(
    "Control",
    "Control",
    "Control",
    "OVA",
    "OVA",
    "OVA",
    "OVA"
  )
)

metadata
summary(asthma[,2:8])
asthma_log <- asthma

asthma_log[,2:8] <- log2(asthma[,2:8] + 1)
summary(asthma_log[,2:8])
boxplot(
  asthma_log[,2:8],
  las = 2,
  main = "Log2(FPKM+1) Distribution"
)
expr_pca <- t(asthma_log[,2:8])
pca <- prcomp(expr_pca, scale. = TRUE)
sum(apply(asthma_log[,2:8], 1, var) == 0)
gene_var <- apply(asthma_log[,2:8], 1, var)

asthma_filtered <- asthma_log[gene_var > 0, ]
dim(asthma_filtered)
expr_pca <- t(asthma_filtered[,2:8])

pca <- prcomp(expr_pca, scale. = TRUE)
plot(
  pca$x[,1],
  pca$x[,2],
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA: Control vs OVA"
)

text(
  pca$x[,1],
  pca$x[,2],
  labels = rownames(expr_pca),
  pos = 3
)
plot(
  pca$x[,1],
  pca$x[,2],
  xlab="PC1",
  ylab="PC2",
  main="PCA: Control vs OVA"
)

text(
  pca$x[,1],
  pca$x[,2],
  labels=rownames(expr_pca),
  pos=3,
  cex=0.8
)
summary(pca)
control_mean <- rowMeans(
  asthma_log[, c("Control72h1","Control72h2","Control72h3")]
)

ova_mean <- rowMeans(
  asthma_log[, c("OVA72h1","OVA72h2","OVA72h3","OVA72h4")]
)
log2FC <- ova_mean - control_mean
results <- data.frame(
  gene_id = asthma_log$gene_id,
  Control_Mean = control_mean,
  OVA_Mean = ova_mean,
  log2FC = log2FC
)
head(results)
upregulated <- results[order(-results$log2FC), ]
head(upregulated, 100)
downregulated <- results[order(results$log2FC), ]
head(downregulated, 20)
top20 <- upregulated[1:20, ]

barplot(
  top20$log2FC,
  names.arg = top20$gene_id,
  las = 2,
  cex.names = 0.7,
  main = "Top Upregulated Genes in OVA Asthma"
)
control_cols <- c("Control72h1","Control72h2","Control72h3")
ova_cols <- c("OVA72h1","OVA72h2","OVA72h3","OVA72h4")
get_pvalue <- function(i) {
  control_vals <- as.numeric(asthma_log[i, control_cols])
  ova_vals <- as.numeric(asthma_log[i, ova_cols])
  
  test <- t.test(ova_vals, control_vals)
  return(test$p.value)
}
pvalues <- sapply(1:nrow(asthma_log), get_pvalue)
results$pvalue <- pvalues
results$padj <- p.adjust(results$pvalue, method = "fdr")
head(results)
