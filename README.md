# U87 LAMP2A RNA-Seq — Starter Repo

This starter contains ready-to-run R scripts that reproduce your GEN3003 pipeline (merge featureCounts → DESeq2 → figures → GO/KEGG).

## Run order (in RStudio Console)
```r
source("scripts/00_merge_featureCounts.R")
source("scripts/01_deseq2_analysis.R")
source("scripts/02_visualizations.R")
source("scripts/03_enrichment.R")
```
> Place your Galaxy featureCounts files (CL1/CL2/CL3) in the project root next to `scripts/` **or** update the file paths inside `00_merge_featureCounts.R` to point to their actual locations.

Outputs are written to `results/` (tables) and `figures/` (plots).