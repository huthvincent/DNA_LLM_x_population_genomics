# ------------------------------------------------------------------------------
# ADNI haplotype association analysis
# Author: Xiaopu Zhou
# Description:
# This script tests the association between each haplotype and AD dementia status
# using logistic regression adjusted for demographic covariates.
# ------------------------------------------------------------------------------


# === Load Required Libraries ===

library(dplyr)

# === Load Data ===

data <- read.table(
  "ADNI_hap_call_1998_with_Pheno_haplotype_Sequence.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

frequency <- read.table(
  "ADNI_hap_call_1087_haplotype_frequency.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# === Prepare Phenotype Data ===

dat_sub <- data %>%
  filter(DX %in% c("Dementia", "CN")) %>%
  mutate(
    DX_bin = as.integer(DX == "Dementia"),
    Month_bl = as.numeric(Month_bl),
    Month = as.numeric(Month),
    PTEDUCAT = as.numeric(PTEDUCAT),
    PTGENDER = as.factor(PTGENDER),
    PTRACCAT = as.factor(PTRACCAT)
  )

# === Helper Function: Fit One Haplotype ===

fit_one_hap <- function(hap) {
  
  df <- dat_sub %>%
    mutate(hap_dose = as.integer(hap1 == hap) + as.integer(hap2 == hap))
  
  if (length(unique(df$hap_dose)) < 2) {
    return(data.frame(
      haplotype = hap,
      beta = NA,
      se = NA,
      OR = NA,
      CI_low = NA,
      CI_high = NA,
      p = NA
    ))
  }
  
  fit <- glm(
    DX_bin ~ hap_dose + Month_bl + Month + PTEDUCAT + PTGENDER + PTRACCAT,
    data = df,
    family = binomial
  )
  
  beta <- coef(fit)["hap_dose"]
  se <- sqrt(vcov(fit)["hap_dose", "hap_dose"])
  OR <- exp(beta)
  CI <- exp(beta + c(-1, 1) * 1.96 * se)
  p <- summary(fit)$coefficients["hap_dose", "Pr(>|z|)"]
  
  data.frame(
    haplotype = hap,
    beta = beta,
    se = se,
    OR = OR,
    CI_low = CI[1],
    CI_high = CI[2],
    p = p
  )
}

# === Run Association Analysis ===

assoc <- bind_rows(lapply(frequency$haplotype, fit_one_hap)) %>%
  mutate(FDR = p.adjust(p, method = "fdr"))

# === Merge with Frequency Table ===

frequency2 <- frequency %>%
  left_join(assoc, by = "haplotype") %>%
  arrange(FDR, p)

# === Write Output ===

write.table(
  frequency2,
  "ADNI_haplotype_frequency_with_GLM_AD_NC.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)