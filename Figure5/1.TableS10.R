# ==============================================================================
# 1.TableS10.R
# Author: Xiaopu Zhou
# Purpose:
#   - Integrate phenotype, covariates, and Evo2 delta scores for UK Biobank
#   - Run logistic regression on AD diagnosis using Evo2 scores
#   - Stratify by ancestry and APOE4 genotype dose
# ==============================================================================

# ==== [0] Set Working Directory ====

# ==== [1] Load Required Data ====

# Load phenotype data (AD vs NC)
data <- read.table("100313_n333791_AllAncestry_unrelated_No.aneuploidy.txt",
                   stringsAsFactors = FALSE, header = TRUE, sep = '\t')
AD <- read.table("AD_cases_n4156.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')
NC <- read.table("NC_cases2_n172850.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')

# Annotate phenotype
data$Pheno <- NA
data$Pheno[data$Participant.ID %in% AD$Participant.ID] <- 'AD'
data$Pheno[data$Participant.ID %in% NC$Participant.ID] <- 'NC'
data <- subset(data, !is.na(Pheno))
data <- data[, c("Participant.ID", "Pheno", "Genetic.ethnic.grouping")]

# Load covariates
anno <- read.table("100094_Age_Sex_clean.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')
anno <- anno[, c("Participant.ID", "Sex", "Age.at.recruitment")]

# Merge phenotype with covariates
data <- merge(data, anno, by = "Participant.ID")

# Load Evo2 score
EVO2_score <- read.table("EVO2_APOE_AD_score_n490155.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')
data <- merge(data, EVO2_score, by.x = "Participant.ID", by.y = "IID")

# Load PCA
PCA <- read.table("100313_PC_n487902.txt", stringsAsFactors = FALSE, header = TRUE)
data <- merge(data, PCA, by = "Participant.ID")

# ==== [2] Prepare Final Dataset ====

# Subset to AD and NC only, adjust Evo2 scores
data_final <- rbind(data[data$Pheno == "AD", ], data[data$Pheno == "NC", ])
data_final$score_leading <- data_final$score_leading - (-0.8497495 * 80000)
data_final$score_alternative <- data_final$score_alternative - (-0.8497495 * 80000)

# ==== [3] Logistic Regression (All Ancestry) ====
test <- glm(
  Pheno == 'AD' ~ score_leading + score_alternative + Sex + Age.at.recruitment +
    PC1 + PC2 + PC3 + PC4 + PC5,
  data = data_final,
  family = binomial
)
summary(test)

