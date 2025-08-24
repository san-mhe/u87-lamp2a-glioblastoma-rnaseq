# scripts/03_enrichment.R

# GO BP + KEGG enrichment and plots (dot/bar) — auto-detects ID type (ENTREZ / ENSEMBL / SYMBOL)

# Install dependencies
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("clusterProfiler","org.Hs.eg.db","AnnotationDbi"))
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
for (pkg in c("dplyr","readr","ggplot2","stringr"))
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
  library(dplyr); library(readr); library(ggplot2); library(stringr)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# Load DEGs (from Top_DEGs_padj0.05_log2FC1.csv)

sig_tbl <- readr::read_csv("results/Top_DEGs_padj0.05_log2FC1.csv", show_col_types = FALSE)

# Choose ID column
candidate_cols <- c("gene","Gene","Geneid","GeneID","ENSEMBL","Ensembl","SYMBOL","Symbol","Row.names","X1")
idcol <- intersect(candidate_cols, names(sig_tbl))
if (length(idcol) == 0) stop("Could not find an ID column in Top_DEGs_padj0.05_log2FC1.csv")
idcol <- idcol[1]

ids_raw <- sig_tbl[[idcol]]

# Detect ID type and normalize for mapping

id_type <- NULL
ids_for_map <- NULL

if (is.numeric(ids_raw)) {
  # Numeric → almost certainly ENTREZ IDs
  id_type <- "ENTREZID"
  ids_for_map <- as.character(unique(ids_raw))
} else {
  ids_chr <- unique(as.character(ids_raw))
  # ENSEMBL with optional version suffix
  if (all(grepl("^ENSG\\d+(\\.\\d+)?$", ids_chr))) {
    id_type <- "ENSEMBL"
    ids_for_map <- sub("\\.\\d+$", "", ids_chr)
  } else {
    # fallback → assume SYMBOL
    id_type <- "SYMBOL"
    ids_for_map <- ids_chr
  }
}

message("Detected ID type: ", id_type)

# Convert to ENTREZ IDs

if (id_type == "ENTREZID") {
  entrez_ids <- unique(ids_for_map)
} else {
  entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys    = ids_for_map,
                                      keytype = id_type,
                                      column  = "ENTREZID",
                                      multiVals = "first") |> unname()
  entrez_ids <- unique(stats::na.omit(entrez_ids))
}

if (length(entrez_ids) == 0) {
  stop("Could not obtain ENTREZ IDs from your input IDs. Check the ID column and type.")
}

# Enrichment analyses
# GO: Biological Process
ego <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_df <- as.data.frame(ego)
readr::write_csv(ego_df, "results/GO_enrichment_results.csv")

# KEGG (human)
ekegg <- enrichKEGG(gene          = entrez_ids,
                    organism      = "hsa",
                    pvalueCutoff  = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekegg_df <- as.data.frame(ekegg)
readr::write_csv(ekegg_df, "results/KEGG_enrichment_results.csv")

# Plots (dot + bar) with robust GeneRatio parsing

parse_ratio <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.character(x)) {
    num <- suppressWarnings(as.numeric(sub("/.*$", "", x)))
    den <- suppressWarnings(as.numeric(sub("^.*/", "", x)))
    return(num / den)
  }
  return(x)
}

go_top  <- ego_df   %>% slice_min(p.adjust, n = 15, with_ties = FALSE)
k_top   <- ekegg_df %>% slice_min(p.adjust, n = 8,  with_ties = FALSE)

# GO plots 
if (nrow(ego_df) > 0) {
  n_top_go <- min(15, nrow(ego_df))
  go_top <- ego_df %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top_go) %>%
    mutate(Description = stringr::str_wrap(Description, width = 55),
           GeneRatioNum = parse_ratio(GeneRatio))
  
  go_bar <- ggplot(go_top, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "red", high = "blue", name = "Adj. p-value") +
    labs(title = "GO Enrichment (Biological Process)", x = "Count", y = NULL) +
    theme_bw(base_size = 14)
  ggsave("figures/GO_BP_barplot.png", go_bar, width = 8, height = 6, dpi = 300)
  
  go_dot <- ggplot(go_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum),
                               size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
    labs(title = "GO Enrichment (Biological Process)", x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 14)
  ggsave("figures/GO_BP_dotplot.png", go_dot, width = 8, height = 6, dpi = 300)
} else {
  message("No GO BP terms passed thresholds; skipping GO plots.")
}


# KEGG plots
if (nrow(ekegg_df) > 0) {
  n_top_kegg <- min(8, nrow(ekegg_df))
  k_top <- ekegg_df %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top_kegg) %>%
    mutate(Description = stringr::str_wrap(Description, width = 55),
           GeneRatioNum = parse_ratio(GeneRatio))
  
  k_bar <- ggplot(k_top, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "red", high = "blue", name = "Adj. p-value") +
    labs(title = "KEGG Pathway Enrichment", x = "Count", y = NULL) +
    theme_bw(base_size = 14)
  ggsave("figures/KEGG_barplot.png", k_bar, width = 8, height = 6, dpi = 300)
  
  k_dot <- ggplot(k_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum),
                             size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
    labs(title = "KEGG Pathway Enrichment", x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 14)
  ggsave("figures/KEGG_dotplot.png", k_dot, width = 8, height = 6, dpi = 300)
} else {
  message("No KEGG pathways passed thresholds; skipping KEGG plots.")
}


# Session info
writeLines(capture.output(sessionInfo()), "results/sessionInfo_03.txt")
