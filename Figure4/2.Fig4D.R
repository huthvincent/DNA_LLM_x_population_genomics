# ==============================================================================
# 2.Fig4D.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merge ADNI haplotype annotations with baseline endophenotype data
#   - Run robust linear regression (lmrob) to test association between Evo2
#     delta scores and cognitive/MRI phenotypes
#   - Analyze both full cohort and subset excluding cognitively normal (CN)
# ==============================================================================

# ==== [0] Set Working Directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Figure3/endophenotype")

# ==== [1] Load and Merge Data ====
hap <- read.table("ADNI_hap_call_1998_with_Pheno.FINAL.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
diag <- read.table("ADNI_baseline_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

hap_anno <- merge(hap, diag, by.x = "ID_update", by.y = "PTID")

# Save merged table
write.table(hap_anno, "ADNI_hap_call_1998_with__baseline_Pheno.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Set reference group for race
hap_anno$PTRACCAT <- relevel(as.factor(hap_anno$PTRACCAT), ref = "White")

# ==== [2] Define APOE4-Negative Subset (optional if needed later) ====
hap_nonE4 <- subset(hap_anno, APOE4_leading_hap == 0 & APOE4_other_hap == 0)

# ==== [3] Load Robust Regression Package ====
library(robustbase)

# ==== [4] Association with All Participants ====

MMSE_test <- lmrob(MMSE ~ delta_leading_hap + delta_other_hap + Month_bl +
                     as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT), data = hap_anno)
summary(MMSE_test)

CDRSB_test <- lmrob(CDRSB ~ delta_leading_hap + delta_other_hap + Month_bl +
                      as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER, data = hap_anno)
summary(CDRSB_test)

ADAS11_test <- lmrob(ADAS11 ~ delta_leading_hap + delta_other_hap + Month_bl +
                       as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT), data = hap_anno)
summary(ADAS11_test)

ADAS13_test <- lmrob(ADAS13 ~ delta_leading_hap + delta_other_hap + Month_bl +
                       as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT), data = hap_anno)
summary(ADAS13_test)

Entorhinal_test <- lmrob(Entorhinal ~ delta_leading_hap + delta_other_hap + Month_bl +
                           as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV, data = hap_anno)
summary(Entorhinal_test)

Hippocampus_test <- lmrob(Hippocampus ~ delta_leading_hap + delta_other_hap + Month_bl +
                            as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV, data = hap_anno)
summary(Hippocampus_test)

WholeBrain_test <- lmrob(WholeBrain ~ delta_leading_hap + delta_other_hap + Month_bl +
                           as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV, data = hap_anno)
summary(WholeBrain_test)

# ==== [5] Association in Patients Only (Exclude CN) ====

MMSE_test <- lmrob(MMSE ~ delta_leading_hap + delta_other_hap + Month_bl +
                     as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
                   data = hap_anno[hap_anno$DX.x != "CN", ])
summary(MMSE_test)

CDRSB_test <- lmrob(CDRSB ~ delta_leading_hap + delta_other_hap + Month_bl +
                      as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER,
                    data = hap_anno[hap_anno$DX.x != "CN", ])
summary(CDRSB_test)

ADAS11_test <- lmrob(ADAS11 ~ delta_leading_hap + delta_other_hap + Month_bl +
                       as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
                     data = hap_anno[hap_anno$DX.x != "CN", ])
summary(ADAS11_test)

ADAS13_test <- lmrob(ADAS13 ~ delta_leading_hap + delta_other_hap + Month_bl +
                       as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT),
                     data = hap_anno[hap_anno$DX.x != "CN", ])
summary(ADAS13_test)

Entorhinal_test <- lmrob(Entorhinal ~ delta_leading_hap + delta_other_hap + Month_bl +
                           as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV,
                         data = hap_anno[hap_anno$DX.x != "CN", ])
summary(Entorhinal_test)

Hippocampus_test <- lmrob(Hippocampus ~ delta_leading_hap + delta_other_hap + Month_bl +
                            as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV,
                          data = hap_anno[hap_anno$DX.x != "CN", ])
summary(Hippocampus_test)

WholeBrain_test <- lmrob(WholeBrain ~ delta_leading_hap + delta_other_hap + Month_bl +
                           as.numeric(Month) + as.numeric(PTEDUCAT) + PTGENDER + as.factor(PTRACCAT) + ICV,
                         data = hap_anno[hap_anno$DX.x != "CN", ])
summary(WholeBrain_test)
