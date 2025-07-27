# ==============================================================================
# 3.Fig4EF.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merge ADNI APOE haplotype data with baseline amyloid PET SUVR measurements
#   - Run linear regression to assess association between Evo2 delta scores
#     and SUVR across full cohort and diagnostic subgroups
# ==============================================================================

# ==== [0] Set Working Directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Figure3/PET_abeta")

# ==== [1] Load and Merge Data ====
hap <- read.table("ADNI_hap_call_1998_with_Pheno.FINAL.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
diag <- read.table("brain_PET_abeta_6mm.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

hap_anno <- merge(hap, diag, by.x = "ID_update", by.y = "PTID")

# Save merged table
write.table(hap_anno, "ADNI_hap_call_1998_with_baseline_PET_abeta_6mm.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Set race reference level
hap_anno$PTRACCAT <- relevel(as.factor(hap_anno$PTRACCAT), ref = "White")

# Define APOE4-negative subset (if needed later)
hap_nonE4 <- subset(hap_anno, APOE4_leading_hap == 0 & APOE4_other_hap == 0)

# ==== [2] Load Robust Regression Package ====
library(robustbase)

# ==== [3] Linear Regression for Whole Cohort ====

SUVR_all <- lm(
  SUMMARY_SUVR ~ delta_leading_hap + delta_other_hap + Month_bl +
    SUMMARY_VOLUME + as.numeric(Month) + PTGENDER,
  data = hap_anno
)
summary(SUVR_all)

# ==== [4] Linear Regression: EMCI & LMCI Subset ====

SUVR_mci <- lm(
  SUMMARY_SUVR ~ delta_leading_hap + delta_other_hap + Month_bl +
    SUMMARY_VOLUME + as.numeric(Month) + PTGENDER,
  data = hap_anno[hap_anno$DX_bl %in% c("EMCI", "LMCI"), ]
)
summary(SUVR_mci)

# ==== [5] Linear Regression: SMC Subset ====

SUVR_smc <- lm(
  SUMMARY_SUVR ~ delta_leading_hap + delta_other_hap + Month_bl +
    SUMMARY_VOLUME + as.numeric(Month) + PTGENDER,
  data = hap_anno[hap_anno$DX_bl == "SMC", ]
)
summary(SUVR_smc)
