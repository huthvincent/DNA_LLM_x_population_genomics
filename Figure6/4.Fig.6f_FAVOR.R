# ------------------------------------------------------------------------------
# APOE/SORL1 FAVOR-Evo2 embedding regression
# Author: Xiaopu Zhou
# Description:
# This script tests associations between individual Evo2 embedding dimensions
# and numeric FAVOR annotation scores, adjusting for genomic region.
# ------------------------------------------------------------------------------


# === Load Required Libraries ===

library(data.table)

# === Load Data ===

data <- fread(
  "APOE_SORL1_387_site_FAVOR_matchFAVOR_allele_Final.txt",
  check.names = FALSE
)

favor_info <- fread(
  "headername.txt",
  check.names = FALSE
)

region <- fread(
  "Evo2_APOE_SORL1_387_site.score.txt",
  check.names = FALSE
)[, .(SNP, Region)]

# === Prepare Headers ===

embedding_header <- names(data)[grepl("^embedding_", names(data))]
favor_header <- favor_info[type == "numeric", names]
covar_header <- "Region"

# === Merge Region Information ===

data <- merge(
  data,
  region,
  by.x = "rsid",
  by.y = "SNP"
)

# === Helper Function ===

run_lm_single <- function(dt, fav, emb, covar_header) {
  
  dt <- dt[, c(fav, emb, covar_header), with = FALSE]
  
  dt[[fav]] <- suppressWarnings(as.numeric(dt[[fav]]))
  dt[[emb]] <- suppressWarnings(as.numeric(dt[[emb]]))
  
  if (sum(is.finite(dt[[emb]])) < 5) return(NULL)
  if (length(unique(dt[[emb]][is.finite(dt[[emb]])])) <= 1) return(NULL)
  
  valid_covars <- character(0)
  
  if (covar_header %in% names(dt)) {
    x <- as.character(dt[[covar_header]])
    x[is.na(x) | x == ""] <- "Missing"
    x <- factor(x)
    
    if (nlevels(x) > 1) {
      dt[[covar_header]] <- x
      valid_covars <- covar_header
    }
  }
  
  complete_cols <- c(fav, emb, valid_covars)
  dt <- dt[complete.cases(dt[, complete_cols, with = FALSE])]
  
  if (nrow(dt) < 10) return(NULL)
  
  form <- as.formula(
    paste(fav, "~", paste(c(emb, valid_covars), collapse = " + "))
  )
  
  fit <- tryCatch(
    lm(form, data = dt),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  coef_res <- summary(fit)$coefficients
  
  if (!(emb %in% rownames(coef_res))) return(NULL)
  
  data.table(
    FAVOR = fav,
    embedding = emb,
    beta = coef_res[emb, "Estimate"],
    SE = coef_res[emb, "Std. Error"],
    t = coef_res[emb, "t value"],
    p = coef_res[emb, "Pr(>|t|)"],
    n = nrow(dt)
  )
}

# === Run Regression Analysis ===

results <- list()
idx <- 1

for (fav in favor_header) {
  
  message("Processing FAVOR: ", fav)
  
  for (emb in embedding_header) {
    
    res <- tryCatch(
      run_lm_single(data, fav, emb, covar_header),
      error = function(e) NULL
    )
    
    if (!is.null(res)) {
      results[[idx]] <- res
      idx <- idx + 1
    }
  }
}

if (length(results) > 0) {
  results <- rbindlist(results, fill = TRUE)
  results[, bonferroni := p.adjust(p, method = "bonferroni")]
  results[, FDR := p.adjust(p, method = "fdr")]
} else {
  results <- data.table()
}

# === Write Output ===

fwrite(
  results,
  "FAVOR_APOE_SORL1_embeding_387_mut_Stage2_single_embedding_lm_results.txt",
  sep = "\t"
)
