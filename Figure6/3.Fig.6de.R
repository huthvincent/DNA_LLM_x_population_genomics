# ------------------------------------------------------------------------------
# FAVOR-Evo2 correlation analysis
# Author: Xiaopu Zhou
# Description:
# This script tests correlations between absolute Evo2 delta scores and numeric
# FAVOR annotation scores.
# ------------------------------------------------------------------------------


# === Load Required Libraries ===

library(readr)
library(dplyr)
library(purrr)

# === Load Data ===

data <- read_tsv(
  "../FAVOR_Evo2_6bnoRC_Inference_result_MATCHallele.2715101.UCSC_SNP.noLCR.txt.gz",
  show_col_types = FALSE
)

header <- read_tsv(
  "headername.txt",
  show_col_types = FALSE
) %>%
  filter(type == "numeric") %>%
  pull(names)

# === Helper Function ===

run_cor <- function(x, y, label) {
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  
  if (n < 3) {
    return(tibble(
      comparison = label,
      n_pairs = n,
      r = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ct <- suppressWarnings(
    cor.test(x, y, method = "pearson")
  )
  
  tibble(
    comparison = label,
    n_pairs = n,
    r = unname(ct$estimate),
    ci_low = ct$conf.int[1],
    ci_high = ct$conf.int[2],
    p_value = ct$p.value
  )
}

# === Run Correlation Analysis ===

x <- abs(data$`7b_noRC_delta_score`)

res_df <- map_dfr(header, function(colname) {
  
  if (!colname %in% names(data)) {
    warning("Column not found: ", colname)
    return(NULL)
  }
  
  y <- suppressWarnings(as.numeric(data[[colname]]))
  
  run_cor(x, y, colname)
}) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# === Write Output ===

write_tsv(
  res_df,
  "FAVOR_Evo2_2715101.UCSC_SNP.noLCR_matchAllele_correlation_results.tsv"
)