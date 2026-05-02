# ------------------------------------------------------------------------------
# Evo2 global annotation analysis
# Author: Xiaopu Zhou
# Description:
# This script compares absolute Evo2 delta scores across functional annotation
# categories using Intergenic variants as the reference group.
# ------------------------------------------------------------------------------

# === Set Working Directory ===

setwd("D:/OneDrive/2026Feb/EVO2_APOE/Revision/7b_global_annotation")

# === Load Required Libraries ===

library(readr)
library(dplyr)

# === Load Data ===

data <- read_tsv(
  "Evo2_6Models_Inference_result.6475639.UCSC_Fullanno.LCR.Jan19.txt.gz",
  show_col_types = FALSE
)

# === Select Columns ===

cols <- c(
  "rsid",
  "BP",
  "Ref",
  "Alt",
  "Alt_Freq",
  "7b_noRC_delta_score",
  "Func.ensGene",
  "Gene.ensGene",
  "GeneDetail.ensGene",
  "category",
  "AAChange.ensGene",
  "phastCons100way",
  "phyloP100way",
  "ExonicFunc.ensGene",
  "repName"
)

data_clean <- data[, cols]

# === Filter Variants ===

data_clean <- data_clean %>%
  filter(
    is.na(repName),
    Ref %in% c("A", "T", "C", "G"),
    Alt %in% c("A", "T", "C", "G")
  ) %>%
  mutate(
    score = abs(`7b_noRC_delta_score`),
    cat = as.factor(category)
  ) %>%
  filter(!is.na(score), !is.na(cat))

# === Summary by Category ===

summary_all <- data_clean %>%
  group_by(cat) %>%
  summarise(
    n = n(),
    mean_score = mean(score),
    median_score = median(score),
    sd_score = sd(score),
    .groups = "drop"
  )

# === Wilcoxon Test Against Intergenic Variants ===

ref <- "Intergenic"
cats_to_test <- setdiff(levels(data_clean$cat), ref)

res <- lapply(cats_to_test, function(g) {
  
  d2 <- data_clean %>%
    filter(cat %in% c(ref, g))
  
  wt <- wilcox.test(score ~ cat, data = d2, exact = FALSE)
  
  data.frame(
    group = g,
    n_ref = sum(d2$cat == ref),
    n_group = sum(d2$cat == g),
    mean_ref = mean(d2$score[d2$cat == ref]),
    mean_group = mean(d2$score[d2$cat == g]),
    mean_diff = mean(d2$score[d2$cat == g]) - mean(d2$score[d2$cat == ref]),
    median_ref = median(d2$score[d2$cat == ref]),
    median_group = median(d2$score[d2$cat == g]),
    median_diff = median(d2$score[d2$cat == g]) - median(d2$score[d2$cat == ref]),
    p_raw = wt$p.value
  )
  
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_raw, method = "fdr")) %>%
  arrange(p_adj)

# === Write Output ===


write.table(
  res,
  "Evo2_global_annotation_wilcoxon_vs_intergenic.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)