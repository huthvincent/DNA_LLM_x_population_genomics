# ==============================================================================
# 8.Fig5D_PheWAS.R
# Author: Xiaopu Zhou
# Purpose:
#   - Performs logistic regression across 1547 phecode traits using UK Biobank
#   - Predictors: score_leading and score_alternative from EVO2 model
#   - Adjusts for covariates: age, sex, top 5 PCs
#   - Outputs effect size, p-value and FDR for each phenotype
# ==============================================================================

# ------------------------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load base phenotype and covariates
# ------------------------------------------------------------------------------
data <- read.table("100313_n333791_AllAncestry_unrelated_No.aneuploidy.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)

anno <- read.table("100094_Age_Sex_clean.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
anno <- anno[, c("Participant.ID", "Sex", "Age.at.recruitment")]

data <- merge(data, anno, by = "Participant.ID")

# ------------------------------------------------------------------------------
# Load EVO2 scores and PCA
# ------------------------------------------------------------------------------
EVO2_score <- read.table("EVO2_APOE_AD_score_n490155.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data <- merge(data, EVO2_score, by.x = "Participant.ID", by.y = "IID")

PCA <- read.table("100313_PC_n487902.txt", header = TRUE, stringsAsFactors = FALSE)
data <- merge(data, PCA, by = "Participant.ID")

# ------------------------------------------------------------------------------
# Load and merge PheWAS data
# ------------------------------------------------------------------------------
library(readr)
phewas <- read_tsv("ICD10_to_Phecode_443702n_1547Pheno.txt")
data <- merge(data, phewas, by.x = "Participant.ID", by.y = "Participant ID", all.x = TRUE)

# Replace missing (NA) values in phecodes with 0
for (column_id in 24:1570) {
  data[[column_id]][is.na(data[[column_id]])] <- 0
}

# ------------------------------------------------------------------------------
# Normalize scores
# ------------------------------------------------------------------------------
data$score_leading <- data$score_leading - (-0.8497495 * 80000)
data$score_alternative <- data$score_alternative - (-0.8497495 * 80000)

data_final <- data

# ------------------------------------------------------------------------------
# Run logistic regression (PheWAS)
# ------------------------------------------------------------------------------
library(robustbase)  # if fallback glmrob is needed
library(RNOmni)

disease_cols <- names(data_final)[24:1570]

results <- data.frame(
  Disease = character(),
  Predictor = character(),
  Beta = numeric(),
  SE = numeric(),
  z_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (disease_id in disease_cols) {
  
  data_final$disease <- data_final[[disease_id]]
  
  formula <- disease ~ score_leading + score_alternative +
    Sex + Age.at.recruitment + PC1 + PC2 + PC3 + PC4 + PC5
  
  try({
    model <- glm(formula, data = data_final, family = binomial, control = glm.control(maxit = 100))
    coef_summary <- summary(model)$coefficients
    
    for (predictor in c("score_leading", "score_alternative")) {
      results <- rbind(results, data.frame(
        Disease = disease_id,
        Predictor = predictor,
        Beta = coef_summary[predictor, "Estimate"],
        SE = coef_summary[predictor, "Std. Error"],
        z_value = coef_summary[predictor, "z value"],
        p_value = coef_summary[predictor, "Pr(>|z|)"]
      ))
    }
  }, silent = TRUE)
}

# ------------------------------------------------------------------------------
# Multiple testing correction
# ------------------------------------------------------------------------------
results <- results[order(results$p_value), ]

results$FDR[results$Predictor == "score_leading"] <- p.adjust(results$p_value[results$Predictor == "score_leading"], method = "fdr")
results$FDR[results$Predictor == "score_alternative"] <- p.adjust(results$p_value[results$Predictor == "score_alternative"], method = "fdr")

# ------------------------------------------------------------------------------
# Output results
# ------------------------------------------------------------------------------
write.table(results, "PheWAS_association_results_1547disease_n331925.raw.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

