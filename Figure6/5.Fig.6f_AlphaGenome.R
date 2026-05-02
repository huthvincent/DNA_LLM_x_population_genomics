# ------------------------------------------------------------------------------
# AlphaGenome-Evo2 embedding regression
# Author: Xiaopu Zhou
# Description:
# This script tests associations between Evo2 embedding dimensions and
# AlphaGenome raw scores across assay types, adjusting for available covariates.
# ------------------------------------------------------------------------------

# === Load Required Libraries ===

library(readr)
library(data.table)

# === Load Data ===

data <- fread(
  "APOE_SORL1_embeding_387_mut_Final.txt",
  check.names = FALSE
)

region <- fread(
  "Evo2_APOE_SORL1_387_site.score.txt",
  check.names = FALSE
)

# === Prepare Data ===

embedding_header <- names(data)[grepl("^embedding_", names(data))]

data <- merge(
  data,
  region,
  by.x = "ID",
  by.y = "SNP"
)

data <- as.data.table(data)

# === Define Covariates and Assay Types ===

covar_header <- c(
  "Region",
  "biosample_name",
  "biosample_type",
  "gene_name",
  "biosample_life_stage",
  "transcription_factor",
  "histone_mark",
  "gtex_tissue"
)

assay_type_list <- c(
  "ATAC",
  "CAGE",
  "CHIP_HISTONE",
  "CHIP_TF",
  "DNASE",
  "PROCAP",
  "RNA_SEQ",
  "SPLICE_JUNCTIONS",
  "SPLICE_SITE_USAGE"
)

# === Helper Function ===

clean_covariates <- function(dt, covars, min_count = 3) {
  
  keep <- character(0)
  
  for (cv in covars) {
    
    if (!cv %in% names(dt)) next
    
    x <- as.character(dt[[cv]])
    x[is.na(x) | x == ""] <- "Missing"
    
    tab <- table(x)
    if (length(tab) <= 1) next
    
    rare_levels <- names(tab)[tab < min_count]
    x[x %in% rare_levels] <- "Other"
    x <- droplevels(factor(x))
    
    if (nlevels(x) <= 1) next
    
    dt[[cv]] <- x
    keep <- c(keep, cv)
  }
  
  list(dt = dt, keep = keep)
}

run_lm_multivariate <- function(dt, y_var, embedding_vars, covars) {
  
  use_cols <- intersect(c(y_var, embedding_vars, covars), names(dt))
  dt <- copy(dt[, use_cols, with = FALSE])
  
  dt[[y_var]] <- suppressWarnings(as.numeric(dt[[y_var]]))
  dt <- dt[is.finite(get(y_var))]
  
  emb_present <- intersect(embedding_vars, names(dt))
  
  for (emb in emb_present) {
    dt[[emb]] <- suppressWarnings(as.numeric(dt[[emb]]))
  }
  
  emb_keep <- emb_present[
    sapply(emb_present, function(emb) {
      x <- dt[[emb]]
      sum(is.finite(x)) >= 5 && length(unique(x[is.finite(x)])) > 1
    })
  ]
  
  if (length(emb_keep) == 0) return(NULL)
  
  dt <- dt[complete.cases(dt[, emb_keep, with = FALSE])]
  if (nrow(dt) < 10) return(NULL)
  
  covar_res <- clean_covariates(dt, covars)
  dt <- covar_res$dt
  valid_covars <- covar_res$keep
  
  fit <- tryCatch(
    lm(
      reformulate(c(emb_keep, valid_covars), response = y_var),
      data = dt
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  coef_res <- summary(fit)$coefficients
  emb_final <- intersect(emb_keep, rownames(coef_res))
  
  if (length(emb_final) == 0) return(NULL)
  
  data.table(
    embedding = emb_final,
    beta = coef_res[emb_final, "Estimate"],
    SE = coef_res[emb_final, "Std. Error"],
    t = coef_res[emb_final, "t value"],
    p = coef_res[emb_final, "Pr(>|t|)"],
    n = nrow(dt),
    n_embedding = length(emb_keep),
    n_covariates = length(valid_covars),
    covariates_used = ifelse(
      length(valid_covars) > 0,
      paste(valid_covars, collapse = ";"),
      NA_character_
    )
  )
}

# === Run Regression Analysis ===

results <- list()
idx <- 1

for (assay_type in assay_type_list) {
  
  message("Processing assay: ", assay_type)
  
  file_in <- paste0("./data_block/AlphaGenome_", assay_type, "_data_block.txt.gz")
  
  if (!file.exists(file_in)) next
  
  alphagenome <- read_tsv(file_in, show_col_types = FALSE)
  alphagenome <- as.data.table(alphagenome)
  
  if (!"variant_id" %in% names(alphagenome)) next
  
  merged_dt <- merge(
    alphagenome,
    data,
    by = "variant_id"
  )
  
  if (!"raw_score" %in% names(merged_dt)) next
  if (nrow(merged_dt) == 0) next
  
  res <- run_lm_multivariate(
    dt = merged_dt,
    y_var = "raw_score",
    embedding_vars = embedding_header,
    covars = covar_header
  )
  
  if (!is.null(res)) {
    res[, assay_type := assay_type]
    results[[idx]] <- res
    idx <- idx + 1
  }
}

# === Adjust P-values ===

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
  "AlphaGenome_APOE_SORL1_embedding_387_mut_raw_score_multivariate_lm_results.txt",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)