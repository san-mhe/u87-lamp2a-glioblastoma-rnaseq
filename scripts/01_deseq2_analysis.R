# scripts/01_deseq2_analysis.R
# KD vs WT with cell-line adjustment; exports DE tables + normalized counts + RDS objects

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2")) if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
for (pkg in c("readr","dplyr","tibble")) if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)

library(DESeq2); library(readr); library(dplyr); library(tibble)

counts_path <- "data/counts/merged_counts_CL1_CL2_CL3.csv"
dir.create("results", showWarnings = FALSE)

#Load counts
counts_raw <- read_csv(counts_path, show_col_types = FALSE)
stopifnot(all(c("Geneid","WT_1","KD_1","WT_2","KD_2","WT_3","KD_3") %in% names(counts_raw)))

counts_mat <- counts_raw %>% column_to_rownames("Geneid") %>% as.matrix()
storage.mode(counts_mat) <- "integer"

#Metadata
coldata <- data.frame(
  row.names = colnames(counts_mat),
  condition = factor(c("WT","KD","WT","KD","WT","KD"), levels = c("WT","KD")),
  cell_line = factor(c("CL1","CL1","CL2","CL2","CL3","CL3"))
)

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ cell_line + condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # low-count filter
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition","KD","WT"))
res_tbl <- as.data.frame(res) %>% rownames_to_column("gene") %>%
  arrange(padj, desc(abs(log2FoldChange)))

#Save
write.csv(as.data.frame(counts_mat), "results/Raw_counts_matrix.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "results/Normalized_counts_all_samples.csv")
write.csv(res_tbl, "results/DESeq2_results_all.csv", row.names = FALSE)

sig_tbl <- res_tbl %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
write.csv(sig_tbl, "results/Top_DEGs_padj0.05_log2FC1.csv", row.names = FALSE)

saveRDS(dds, "results/dds.rds")
saveRDS(res, "results/res.rds")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_01.txt")
