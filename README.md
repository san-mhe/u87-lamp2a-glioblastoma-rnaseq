U87 Glioblastoma RNA-Seq (LAMP2A Knockdown vs WT)

Differential Expression and Pathway Enrichment Analysis using DESeq2

Project Overview

This project analyzes RNA-Seq data from U87 glioblastoma cells, comparing wild-type (WT) against LAMP2A knockdown (KD) samples.

Why LAMP2A?
LAMP2A is a key receptor in chaperone-mediated autophagy (CMA). CMA helps glioblastoma cells adapt to stress and evade apoptosis. Knockdown of LAMP2A is expected to disrupt survival pathways and reveal novel therapeutic vulnerabilities.

Goals:

Identify differentially expressed genes (DEGs) between KD and WT cells.

Perform GO (Biological Process) and KEGG pathway enrichment to understand functional changes.

Visualize results with PCA, volcano plots, heatmaps, and enrichment plots.

Data Retrieval

RNA-Seq count data was retrieved from NCBI SRA (Galaxy-processed featureCounts tables):

Samples used:

WT replicates (3 cell lines)

KD replicates (3 cell lines)

Source: Chaperone-mediated autophagy in glioblastoma, SRA BioProject (example: SRR29525293, SRR29525290
)

Raw reads were aligned and quantified using Galaxy featureCounts → imported into R.

Pipeline

The analysis is structured into 3 R scripts:

01_deseq2.R – Import counts, build DESeq2 dataset, normalize, run differential expression.

Outputs: DESeq2_results_all.csv, Top_DEGs_padj0.05_log2FC1.csv

02_visualization.R – Generate QC & DEG visualizations:

PCA plot

MA plot

Volcano plot

Heatmap of top 40 DEGs

03_enrichment.R – Functional analysis of DEGs:

GO (Biological Process) enrichment

KEGG pathway enrichment

Dotplots and barplots

How to run:

# Clone repo
git clone https://github.com/san-mhe/u87-lamp2a-glioblastoma-rnaseq.git
cd u87-lamp2a-glioblastoma-rnaseq

# Open R / RStudio, then run scripts in order:
source("scripts/01_deseq2.R")
source("scripts/02_visualization.R")
source("scripts/03_enrichment.R")

Results

PCA & QC

WT and KD samples separate clearly, confirming biological differences.

<img src="figures/PCA.png" width="450"/> <img src="figures/MA_plot.png" width="450"/>

Differential Expression

~2,000 DEGs detected (padj < 0.05).

Top upregulated genes in KD include CXCL8, ICAM1, linked to immune cell recruitment.

Top downregulated genes include MEST, HGF, tied to tumor growth and signaling.

<img src="figures/Volcano.png" width="500"/>

GO Biological Processes

Enriched processes include:

Leukocyte migration & chemotaxis

Cell adhesion & junction assembly

Synaptic signaling regulation

These suggest immune microenvironment remodeling upon LAMP2A knockdown.

<img src="figures/GO_BP_barplot.png" width="500"/>

KEGG Pathways

Enriched pathways include:

Cytoskeleton in muscle cells

Dilated/Hypertrophic cardiomyopathy

Cell adhesion molecules

Focal adhesion

This reflects cytoskeletal reorganization and adhesion signaling disruption.

<img src="figures/KEGG_barplot.png" width="500"/>

Biological Interpretation

CMA inhibition (via LAMP2A KD) shifts glioblastoma cells toward immune-related and adhesion pathways, potentially making them more vulnerable to immune clearance.

Targeting CMA in glioblastoma could synergize with immunotherapy or adhesion-targeting drugs.

Dependencies

R ≥ 4.3

Packages: DESeq2, pheatmap, ggplot2, dplyr, clusterProfiler, org.Hs.eg.db, enrichplot

Install with:

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db","enrichplot"))
install.packages(c("ggplot2","dplyr","pheatmap","readr","stringr"))
