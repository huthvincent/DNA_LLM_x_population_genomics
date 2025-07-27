# ==============================================================================
# 4.Fig5A_SWI.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merges pangenome delta scores with MRI susceptibility metrics
#   - Filters and normalizes SWI regions
#   - Performs robust linear regression (lmrob) to assess EVO2 score effects
#   - Stratifies results by predictor and MRI modality (T2* vs susceptibility)
#   - Adjusts p-values using FDR and outputs result table
# ==============================================================================

# Set working directory
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure5_UKB/03.SWI")

# ------------------------------------------------------------------------------
# Load and merge phenotype, covariate, and MRI data
# ------------------------------------------------------------------------------

# Load core sample metadata (unrelated, no aneuploidy)
data <- read.table("100313_n333791_AllAncestry_unrelated_No.aneuploidy.txt",
                   stringsAsFactors = FALSE, header = TRUE, sep = '\t')

# Load Age and Sex annotations
anno <- read.table("100094_Age_Sex_clean.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')
anno <- anno[, c("Participant.ID", "Sex", "Age.at.recruitment")]
data <- merge(data, anno, by = "Participant.ID")

# Load EVO2 scores
EVO2_score <- read.table("EVO2_APOE_AD_score_n490155.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')
data <- merge(data, EVO2_score, by.x = "Participant.ID", by.y = "IID")

# Load genetic PCs
PCA <- read.table("100313_PC_n487902.txt", stringsAsFactors = FALSE, header = TRUE)
data <- merge(data, PCA, by = "Participant.ID")

# Load SWI MRI data (T2* and susceptibility)
SWI <- read.csv("109.csv", stringsAsFactors = FALSE, header = TRUE)

# Filter for valid T2* measures (right nucleus accumbens)
SWI <- SWI[!is.na(SWI$Median.T2star.in.accumbens..right....Instance.2), ]

# Optional: total gray volume (if needed for QC)
SWI$Sum_Gray <- apply(SWI[, grepl("Volume", names(SWI))], 1, sum)

# Merge MRI data into main dataset
data <- merge(data, SWI, by = "Participant.ID")

# Keep selected columns (skip derived Volume columns)
data <- data[, c(1, 10:23, 25:56)]

# Optionally restrict to older participants (commented)
# data_final <- data[data$Age.at.recruitment > 60, ]
data_final <- data

# ------------------------------------------------------------------------------
# Normalize EVO2 score magnitudes
# ------------------------------------------------------------------------------

# Offset to center the scores around 0 (same as gray matter analysis)
data_final$score_leading <- data_final$score_leading - (-0.8497495 * 80000)
data_final$score_alternative <- data_final$score_alternative - (-0.8497495 * 80000)

# ------------------------------------------------------------------------------
# Association Testing: Robust regression using lmrob
# ------------------------------------------------------------------------------

library(robustbase)  # robust linear model
library(RNOmni)      # Rank Normalization

# Identify SWI metric columns (manually or via index)
region_cols <- names(data_final)[16:47]

# Storage for result summary
results <- data.frame(
  Region = character(),
  Predictor = character(),
  Beta = numeric(),
  SE = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through SWI brain metrics
for (region in region_cols) {
  # Remove NA rows
  data_tmp <- data_final[!is.na(data_final[[region]]), ]
  
  # Rank normal transformation of MRI region
  data_tmp$Region_norm <- RankNorm(data_tmp[[region]])
  
  # Linear model formula
  formula <- Region_norm ~ score_leading + score_alternative +
    Sex + Age.at.recruitment + PC1 + PC2 + PC3 + PC4 + PC5
  
  # Fit robust model with try-catch
  try({
    model <- lmrob(formula, data = data_tmp, k.max = 2000)
    coef_summary <- summary(model)$coefficients
    
    # Append both scores' results
    for (predictor in c("score_leading", "score_alternative")) {
      results <- rbind(results, data.frame(
        Region = region,
        Predictor = predictor,
        Beta = coef_summary[predictor, "Estimate"],
        SE = coef_summary[predictor, "Std. Error"],
        t_value = coef_summary[predictor, "t value"],
        p_value = coef_summary[predictor, "Pr(>|t|)"]
      ))
    }
  }, silent = TRUE)
}

# ------------------------------------------------------------------------------
# Adjust p-values using FDR separately for T2* and susceptibility metrics
# ------------------------------------------------------------------------------

# Order by p-value
results <- results[order(results$p_value), ]

# FDR correction for T2* metrics
results$FDR[results$Predictor == "score_leading" & grepl("T2star", results$Region)] <-
  p.adjust(results$p_value[results$Predictor == "score_leading" & grepl("T2star", results$Region)], method = "fdr")
results$FDR[results$Predictor == "score_alternative" & grepl("T2star", results$Region)] <-
  p.adjust(results$p_value[results$Predictor == "score_alternative" & grepl("T2star", results$Region)], method = "fdr")

# FDR correction for susceptibility metrics
results$FDR[results$Predictor == "score_leading" & grepl("susceptibility", results$Region)] <-
  p.adjust(results$p_value[results$Predictor == "score_leading" & grepl("susceptibility", results$Region)], method = "fdr")
results$FDR[results$Predictor == "score_alternative" & grepl("susceptibility", results$Region)] <-
  p.adjust(results$p_value[results$Predictor == "score_alternative" & grepl("susceptibility", results$Region)], method = "fdr")

# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------

# Preview top associations
head(results)

# Write full result table
write.table(results, "SWI_gray_volume_association_results_27749.raw.txt", quote = FALSE, row.names = FALSE, sep = '\t')
