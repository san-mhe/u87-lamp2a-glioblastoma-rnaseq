# scripts/02_visualizations.R

#PCA, MA, Volcano, Heatmap (top 40) mirroring your originals
for (pkg in c("DESeq2","ggplot2","pheatmap","dplyr","readr","RColorBrewer","tibble"))
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)

library(DESeq2); library(ggplot2); library(pheatmap); library(dplyr); library(readr)
library(RColorBrewer); library(tibble)

dir.create("figures", showWarnings = FALSE)

dds <- readRDS("results/dds.rds")
res <- readRDS("results/res.rds")
sig_tbl <- read_csv("results/Top_DEGs_padj0.05_log2FC1.csv", show_col_types = FALSE)

#PCA (VST)
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
pca_df <- cbind(as.data.frame(pca$x[,1:2]), as.data.frame(colData(vsd))) %>% rownames_to_column("sample")
p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = cell_line)) +
  geom_point(size = 4) +
  labs(title = "PCA on VST counts", x = paste0("PC1 (",pct[1],"%)"), y = paste0("PC2 (",pct[2],"%)"))
ggsave("figures/PCA.png", p, width = 7, height = 6, dpi = 300)

#MA
png("figures/MA Plot.png", width = 1600, height = 1200, res = 200)
plotMA(res, ylim = c(-5,5), main="DESeq2 MA: KD vs WT (cell line adjusted)")
abline(h = 0, col = "grey")
dev.off()

#Volcano
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>%
  mutate(status = case_when(
    !is.na(padj) & padj < 0.05 & log2FoldChange >= 1  ~ "Up",
    !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
    TRUE ~ "NotSig"))
v <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = status)) +
  geom_point(alpha=.7, size=1.6) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  scale_color_manual(values=c("Down"="#377eb8","Up"="#e41a1c","NotSig"="grey70")) +
  labs(title="Volcano: KD vs WT (cell line adjusted)", x="log2FC", y="-log10 FDR")
ggsave("figures/Volcano Plot.png", v, width = 7.5, height = 6.5, dpi = 300)

#Heatmap (Top 40 by FDR then |log2FC|)

sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ]
sig <- sig[order(sig$padj), ]
top_genes <- rownames(sig)[1:min(40, nrow(sig))]

mat <- assay(vsd)[top_genes, , drop = FALSE]

cd <- as.data.frame(colData(vsd))
ord <- order(cd$condition, cd$cell_line)
mat <- mat[, ord, drop=FALSE]

ann_col <- data.frame(Condition = cd$condition[ord], CellLine = cd$cell_line[ord])
rownames(ann_col) <- colnames(mat)

palette <- colorRampPalette(rev(brewer.pal(9,"RdBu")))(255)
pheatmap(mat,
         color = palette, cluster_cols = FALSE, annotation_col = ann_col,
         main = "Top 40 Differentially Expressed Genes",
         filename = "figures/Heatmap_top_DEGs.png",
         width = 7, height = 9)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_02.txt")
