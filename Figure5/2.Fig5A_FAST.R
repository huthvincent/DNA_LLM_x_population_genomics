# ==============================================================================
# 2.Fig5A_FAST.R
# Author: Xiaopu Zhou
# Purpose:
#   - Evaluate association between Evo2 delta scores and MRI-based gray matter volume
#   - Use robust regression across multiple brain regions from UK Biobank
#   - Normalize brain region volumes and control for global gray matter
#   - Output coefficient estimates and FDR-adjusted p-values
# ==============================================================================

# ==== [0] Set Working Directory ====

# ==== [1] Load Required Libraries ====
library(robustbase)  # for lmrob
library(RNOmni)       # for RankNorm

# ==== [2] Load Participant and Score Data ====
data <- read.table("100313_n333791_AllAncestry_unrelated_No.aneuploidy.txt", header = TRUE, sep = '\t')
anno <- read.table("100094_Age_Sex_clean.txt", header = TRUE, sep = '\t')[, c("Participant.ID", "Sex", "Age.at.recruitment")]
data <- merge(data, anno, by = "Participant.ID")

EVO2_score <- read.table("EVO2_APOE_AD_score_n490155.txt", header = TRUE, sep = '\t')
data <- merge(data, EVO2_score, by.x = "Participant.ID", by.y = "IID")

PCA <- read.table("100313_PC_n487902.txt", header = TRUE)
data <- merge(data, PCA, by = "Participant.ID")

# ==== [3] Load and Clean MRI Volume Data ====
MRI <- read.csv("1101.csv", header = TRUE)
MRI <- MRI[!is.na(MRI$Volume.of.grey.matter.in.X.Cerebellum..right....Instance.2), ]

# Sum gray matter volumes and filter out outliers
MRI$Sum_Gray <- apply(MRI[, grepl("Volume", names(MRI))], 1, sum)
lower_q <- quantile(MRI$Sum_Gray, probs = 0.01, na.rm = TRUE)
upper_q <- quantile(MRI$Sum_Gray, probs = 0.99, na.rm = TRUE)
MRI_filtered <- MRI[MRI$Sum_Gray >= lower_q & MRI$Sum_Gray <= upper_q, ]

cat("Removed", nrow(MRI) - nrow(MRI_filtered), "outliers\n")

# ==== [4] Merge MRI and Finalize Dataset ====
data <- merge(data, MRI_filtered, by = 'Participant.ID')
data <- data[, c(1, 10:163)]  # Drop unnecessary original columns

# Center scores for stability
data$score_leading <- data$score_leading - (-0.8497495 * 80000)
data$score_alternative <- data$score_alternative - (-0.8497495 * 80000)

data_final <- data

# ==== [5] Run Region-wise Robust Regression ====

# Identify brain region volume columns
region_cols <- grep("Volume", names(data_final), value = TRUE)

results <- data.frame(
  Region = character(),
  Predictor = character(),
  Beta = numeric(),
  SE = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each region
for (region in region_cols) {
  data_tmp <- data_final
  
  # Rank-normalize both target region and global volume
  data_tmp$Region_norm <- RankNorm(data_tmp[[region]])
  data_tmp$Sum_Gray_norm <- RankNorm(data_tmp$Sum_Gray)
  
  # Build model formula
  formula <- Region_norm ~ score_leading + score_alternative +
    Sex + Age.at.recruitment + PC1 + PC2 + PC3 + PC4 + PC5 + Sum_Gray_norm
  
  try({
    model <- lmrob(formula, data = data_tmp)
    coef_summary <- summary(model)$coefficients
    
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

# ==== [6] Adjust for Multiple Testing (FDR) ====
results <- results[order(results$p_value), ]
results$FDR[results$Predictor == "score_leading"] <- p.adjust(results$p_value[results$Predictor == "score_leading"], method = 'fdr')
results$FDR[results$Predictor == "score_alternative"] <- p.adjust(results$p_value[results$Predictor == "score_alternative"], method = 'fdr')

# ==== [7] Output Results ====
write.table(results, "MRI_gray_volume_association_results_29783.txt", quote = FALSE, row.names = FALSE, sep = '\t')
