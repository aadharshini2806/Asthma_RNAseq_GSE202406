# Calculate -log10 adjusted p-value
results$negLogP <- -log10(results$padj)

# Categorize genes
results$Significance <- "Not Significant"

results$Significance[
  results$padj < 0.05 & results$log2FC > 1
] <- "Upregulated"

results$Significance[
  results$padj < 0.05 & results$log2FC < -1
] <- "Downregulated"

# Volcano plot
plot(
  results$log2FC,
  results$negLogP,
  pch = 20,
  col = ifelse(
    results$Significance == "Upregulated", "red",
    ifelse(results$Significance == "Downregulated", "blue", "grey")
  ),
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  main = "Volcano Plot"
)

abline(v = c(-1, 1), col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)
sig_genes <- subset(
  results,
  padj < 0.05 & abs(log2FC) > 1
)

nrow(sig_genes)
up_genes <- sig_genes[
  order(-sig_genes$log2FC),
]

down_genes <- sig_genes[
  order(sig_genes$log2FC),
]
top20 <- head(up_genes, 20)

top20
heatmap_data <- asthma_log[
  asthma_log$gene_id %in% top20$gene_id,
  c(control_cols, ova_cols)
]

rownames(heatmap_data) <- asthma_log$gene_id[
  asthma_log$gene_id %in% top20$gene_id
]

heatmap_matrix <- as.matrix(heatmap_data)

heatmap(
  heatmap_matrix,
  scale = "row",
  main = "Top 20 Upregulated Genes"
)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library(clusterProfiler)
library(org.Mm.eg.db)
genes <- sig_genes$gene_id

gene_map <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

head(gene_map)
dim(gene_map)
go_bp <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(go_bp))
dotplot(go_bp, showCategory = 15)

go_mf <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(go_mf))
dotplot(go_mf, showCategory = 15)
go_cc <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(as.data.frame(go_cc))
dotplot(go_cc, showCategory = 15)
write.csv(
  as.data.frame(go_bp),
  "results/GO_BP_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(go_mf),
  "results/GO_MF_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(go_cc),
  "results/GO_CC_results.csv",
  row.names = FALSE
)