# ------------------------------------------------------------------------------
# ADNI haplotype Evo2 correlation analysis
# Author: Xiaopu Zhou
# Description:
# This script tests the correlation between haplotype AD effect sizes and Evo2
# haplotype scores using Spearman correlation and bootstrap resampling.
# ------------------------------------------------------------------------------

# === Load Required Libraries ===

library(dplyr)

# === Load Data ===

beta_AD <- read.table(
  "ADNI_haplotype_frequency_with_GLM_AD_NC.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

EVO2 <- read.table(
  "ADNI_Haplotype_score_match.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# === Merge Data ===

beta_AD <- merge(beta_AD, EVO2, by = "haplotype")
beta_AD <- na.omit(beta_AD)

beta_AD_selected <- beta_AD %>%
  filter(CN_freq > 0.007)

# === Spearman Correlation ===

cor.test(
  as.numeric(beta_AD_selected$beta),
  as.numeric(beta_AD_selected$expected_delta),
  method = "spearman"
)

# === Bootstrap Correlation Test ===

set.seed(1)

B <- 500
prop <- 0.8
method <- "spearman"

x <- as.numeric(beta_AD_selected$beta)
y <- as.numeric(beta_AD_selected$expected_delta)

ok <- is.finite(x) & is.finite(y)
x <- x[ok]
y <- y[ok]

n <- length(x)
m <- floor(n * prop)

rho_boot <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample.int(n, m, replace = FALSE)
  rho_boot[b] <- cor(x[idx], y[idx], method = method)
}

CI <- quantile(rho_boot, c(0.025, 0.975), na.rm = TRUE)

p_boot <- 2 * min(
  mean(rho_boot <= 0, na.rm = TRUE),
  mean(rho_boot >= 0, na.rm = TRUE)
)

rho_full <- cor(x, y, method = method)

# === Output Results ===

c(
  rho_full = rho_full,
  CI_low = CI[1],
  CI_high = CI[2],
  p_boot = p_boot
)