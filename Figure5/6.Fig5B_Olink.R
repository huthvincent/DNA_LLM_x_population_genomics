# ==============================================================================
# 6.Fig5B_Olink.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merges Olink proteomics data with sample metadata
#   - Applies robust linear regression to assess the association between
#     EVO2 scores and protein levels, adjusting for covariates
#   - Performs FDR correction and outputs summary statistics
# ==============================================================================

# Set working directory
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure5_UKB/04.Olink")

# Load libraries
library(robustbase)   # for lmrob
library(RNOmni)       # for RankNorm
library(readr)        # for read_tsv

# ------------------------------------------------------------------------------
# Load core data
# ------------------------------------------------------------------------------

# Sample metadata
meta <- read.table("100313_n333791_AllAncestry_unrelated_No.aneuploidy.txt",
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Age and sex covariates
anno <- read.table("100094_Age_Sex_clean.txt",
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anno <- anno[, c("Participant.ID", "Sex", "Age.at.recruitment")]

# Merge metadata and covariates
data <- merge(meta, anno, by = "Participant.ID")

# EVO2 scores
EVO2_score <- read.table("EVO2_APOE_AD_score_n490155.txt",
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data <- merge(data, EVO2_score, by.x = "Participant.ID", by.y = "IID")

# PCA covariates
PCA <- read.table("100313_PC_n487902.txt", header = TRUE, stringsAsFactors = FALSE)
data <- merge(data, PCA, by = "Participant.ID")

# ------------------------------------------------------------------------------
# Load Olink proteomics data and merge
# ------------------------------------------------------------------------------

Olink <- read_tsv("olink_clean_53012n_2922protein.txt")
data <- merge(data, Olink, by.x = "Participant.ID", by.y = "eid")

# Keep necessary columns (exclude duplicated participant columns if any)
data <- data[, c(1, 10:ncol(data))]

# Optional filter by age if needed
# data_final <- data[data$Age.at.recruitment > 60, ]

data_final <- data

# Adjust scores to remove offset
score_offset <- -0.8497495 * 80000
data_final$score_leading     <- data_final$score_leading     - score_offset
data_final$score_alternative <- data_final$score_alternative - score_offset

# ------------------------------------------------------------------------------
# Association analysis
# ------------------------------------------------------------------------------

# Get Olink protein names
region_cols <- names(data_final)[24:ncol(data_final)]

# Prepare results container
results <- data.frame(
  Region = character(),
  Predictor = character(),
  Beta = numeric(),
  SE = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each protein
for (region in region_cols) {
  
  data_tmp <- data_final[!is.na(data_final[[region]]), ]
  
  # Apply Rank Normalization
  data_tmp$Region_norm <- RankNorm(data_tmp[[region]])
  
  # Linear model formula
  formula <- Region_norm ~ score_leading + score_alternative +
    Sex + Age.at.recruitment + PC1 + PC2 + PC3 + PC4 + PC5
  
  # Robust regression with fallback
  try({
    model <- lmrob(formula, data = data_tmp, k.max = 2000, maxit.scale = 2000)
    coef_summary <- summary(model)$coefficients
    
    for (predictor in c("score_leading", "score_alternative")) {
      results <- rbind(results, data.frame(
        Region    = region,
        Predictor = predictor,
        Beta      = coef_summary[predictor, "Estimate"],
        SE        = coef_summary[predictor, "Std. Error"],
        t_value   = coef_summary[predictor, "t value"],
        p_value   = coef_summary[predictor, "Pr(>|t|)"]
      ))
    }
  }, silent = TRUE)
}

# FDR correction
results <- results[order(results$p_value), ]
results$FDR <- NA
results$FDR[results$Predictor == "score_leading"]     <- p.adjust(results$p_value[results$Predictor == "score_leading"], method = "fdr")
results$FDR[results$Predictor == "score_alternative"] <- p.adjust(results$p_value[results$Predictor == "score_alternative"], method = "fdr")

# Output top results preview
head(results)

# Save results
write.table(results, "Olink_association_results_36021.raw.txt", quote = FALSE, row.names = FALSE, sep = "\t")

