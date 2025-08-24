# scripts/00_merge_featureCounts.R
# Merge featureCounts outputs from GALAXY into a 6-column raw counts matrix.

if (!requireNamespace("readr", quietly=TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
library(readr); library(dplyr)


fn_CL1_WT <- "Galaxy120_featureCounts_on_data_42_Counts_WT_U187_1.tabular"
fn_CL1_KD <- "Galaxy122_featureCounts_on_data_43_Counts_KD_U187_1.tabular"
fn_CL2_combined <- "Galaxy46-[Cut on data 45].tabular"  # first 3 cols: Geneid, WT_2, KD_2
fn_CL3_WT <- "Galaxy187_featureCounts_on_data_18_Counts_WT_UT187_3.tabular"
fn_CL3_KD <- "Galaxy188_featureCounts_on_data_21_Counts_KD_UT187_3.tabular"

# Standard featureCounts: keep Geneid + last numeric column (counts)
read_fc_lastcol <- function(file) {
  x <- read_tsv(file, col_types = cols(.default = col_double(), Geneid = col_character()))
  tibble(Geneid = x$Geneid, counts = x[[ncol(x)]])
}

# CL2 combined file: first three columns are Geneid, WT, KD
read_cl2_combo <- function(file, wt_name = "WT_2", kd_name = "KD_2") {
  x <- read_tsv(file, col_types = cols(.default = col_double(), Geneid = col_character()))
  colnames(x)[1:3] <- c("Geneid", wt_name, kd_name)
  x[, 1:3]
}

# Read each dataset
CL1_WT <- read_fc_lastcol(fn_CL1_WT) %>% rename(WT_1 = counts)
CL1_KD <- read_fc_lastcol(fn_CL1_KD) %>% rename(KD_1 = counts)
CL2    <- read_cl2_combo(fn_CL2_combined, wt_name = "WT_2", kd_name = "KD_2")
CL3_WT <- read_fc_lastcol(fn_CL3_WT) %>% rename(WT_3 = counts)
CL3_KD <- read_fc_lastcol(fn_CL3_KD) %>% rename(KD_3 = counts)

# Merge by Geneid, fill NAs, collapse duplicates
merged <- CL1_WT %>%
  full_join(CL1_KD, by = "Geneid") %>%
  full_join(CL2,    by = "Geneid") %>%
  full_join(CL3_WT, by = "Geneid") %>%
  full_join(CL3_KD, by = "Geneid")

merged[is.na(merged)] <- 0

merged <- merged %>%
  group_by(Geneid) %>%
  summarise(across(everything(), sum), .groups = "drop")

dir.create("data/counts", recursive = TRUE, showWarnings = FALSE)
write.csv(merged, "data/counts/merged_counts_CL1_CL2_CL3.csv", row.names = FALSE)

message("Saved: data/counts/merged_counts_CL1_CL2_CL3.csv")
