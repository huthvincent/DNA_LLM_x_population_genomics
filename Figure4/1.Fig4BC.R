# ==============================================================================
# 1.Fig4BC.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merge ADNI haplotype annotations with latest diagnosis and metadata
#   - Generate demographic summaries across diagnosis groups
#   - Run logistic regression to test association between Evo2 delta scores and
#     Alzheimer's disease diagnosis (Dementia vs CN / MCI), stratified by APOE4 dose
# ==============================================================================

# ==== [0] Set Working Directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Figure3/AD_diagonosis")

# ==== [1] Load Data ====
hap <- read.table("ADNI_hap_call_final_updateAPOE_reORDER.txt", header = TRUE, stringsAsFactors = FALSE)
diag <- read.table("ADNI_latest_diag.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge genotype + phenotype data
hap_anno <- merge(hap, diag, by.x = "ID_update", by.y = "PTID")

# Save merged table
write.table(hap_anno, "ADNI_hap_call_1998_with_Pheno.FINAL.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Relevel race category to set "White" as baseline
hap_anno$PTRACCAT <- relevel(as.factor(hap_anno$PTRACCAT), ref = "White")

# ==== [2] Demographic Summary ====

# Diagnosis count
table(hap_anno$DX)

# Gender × Diagnosis
table(hap_anno$DX, hap_anno$PTGENDER)

# Age at visit
hap_anno$AGE_update <- hap_anno$AGE + hap_anno$Years_bl
aggregate(AGE_update ~ DX, data = hap_anno, FUN = mean)
aggregate(AGE_update ~ DX, data = hap_anno, FUN = sd)

# Education
aggregate(PTEDUCAT ~ DX, data = hap_anno, FUN = mean)
aggregate(PTEDUCAT ~ DX, data = hap_anno, FUN = sd)

# APOE allele carriers × Diagnosis
table((hap_anno$APOE4_leading_hap + hap_anno$APOE4_other_hap) > 0, hap_anno$DX)
table((hap_anno$APOE2_leading_hap + hap_anno$APOE2_other_hap) > 0, hap_anno$DX)

# ==== [3] Logistic Regression: Dementia vs CN ====
library(robustbase)

# Full model
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "CN")),
  family = binomial
)
summary(test_AD)

# APOE4 = 0 (non-carriers)
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "CN") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 0),
  family = binomial
)
summary(test_AD)

# APOE4 = 1 (heterozygous)
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap * delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "CN") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 1),
  family = binomial
)
summary(test_AD)

# APOE4 = 2 (homozygous)
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "CN") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 2),
  family = binomial
)
summary(test_AD)

# Binary haplotype effect among APOE4-negative
test_AD <- glm(
  DX == "Dementia" ~ as.numeric(delta_leading_hap > 0) + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "CN") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 0),
  family = binomial
)
summary(test_AD)

# ==== [4] Logistic Regression: Dementia vs MCI ====

# Full model
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "MCI")),
  family = binomial
)
summary(test_AD)

# APOE4 = 0
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "MCI") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 0),
  family = binomial
)
summary(test_AD)

# APOE4 = 1
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap * delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "MCI") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 1),
  family = binomial
)
summary(test_AD)

# APOE4 = 2
test_AD <- glm(
  DX == "Dementia" ~ delta_leading_hap + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "MCI") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 2),
  family = binomial
)
summary(test_AD)

# Binary leading haplotype effect in APOE4=0 group
test_AD <- glm(
  DX == "Dementia" ~ as.numeric(delta_leading_hap > 0) + delta_other_hap + Month_bl +
    as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
  data = subset(hap_anno, DX %in% c("Dementia", "MCI") & 
                  (APOE4_leading_hap + APOE4_other_hap) == 0),
  family = binomial
)
summary(test_AD)
