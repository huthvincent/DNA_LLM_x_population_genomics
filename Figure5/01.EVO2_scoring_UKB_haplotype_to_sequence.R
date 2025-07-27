# ==============================================================================
# 01.EVO2_scoring_UKB_haplotype_to_sequence.R
# Author: Xiaopu Zhou
# Purpose:
#   - Generate per-haplotype FASTA sequences for APOE region (UK Biobank WGS)
#   - Input: Clean reference, haplotype-resolved SNP table, haplotype dictionary
#   - Output: One FASTA per unique haplotype (e.g., 490558_Hap1)
# ==============================================================================

# ==== [0] Set Working Directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Analysis/Task4")

# ==== [1] Load Required Libraries ====
library(readr)

# ==== [2] Load Reference Sequence ====
ref <- read_tsv("../APOE_ref_44841743_44921743_GRCh38_clean.txt", col_names = FALSE)$X1
cat("Reference sequence length:", nchar(ref), "\n")  # should equal 80001

# ==== [3] Load Variant Table and Haplotype Dictionary ====
snp_table <- read_tsv("./All_site_genotype/UKB_490558_203_haplotype_table.txt")
hap_table <- read_tsv("UKB_AD203_haplotype.APOE.dict")  # contains: ID, hap_type

# Validate column match
sample_cols <- names(snp_table)[6:ncol(snp_table)]
if (!all(sample_cols == hap_table$ID)) {
  warning("Mismatch between SNP sample names and haplotype dictionary.")
}

# Replace sample IDs with haplotype labels
snp_table_unique <- snp_table
colnames(snp_table_unique)[6:ncol(snp_table_unique)] <- hap_table$hap_type
snp_table_unique <- snp_table_unique[, !duplicated(colnames(snp_table_unique))]

# ==== [4] Prepare to Generate FASTA Sequences ====
ref_seq <- ref
haplotypes <- colnames(snp_table_unique)[6:ncol(snp_table_unique)]
ref_start <- 44841743

# Output directory
dir.create("fasta_output", showWarnings = FALSE)

# ==== [5] Loop Over Haplotype Columns ====
for (hap in haplotypes) {
  hap_seq <- ref_seq  # Initialize with reference sequence
  
  for (i in 1:nrow(snp_table_unique)) {
    pos <- snp_table_unique$POS[i]
    ref_base <- snp_table_unique$REF[i]
    alt_bases <- unlist(strsplit(snp_table_unique$ALT[i], ","))  # multiallelic support
    allele_idx <- as.integer(snp_table_unique[i, hap])
    
    idx <- pos - ref_start + 1  # position in ref string
    ref_len <- nchar(ref_base)
    
    # Check reference base match
    ref_in_seq <- substring(hap_seq, idx, idx + ref_len - 1)
    if (ref_in_seq != ref_base) {
      warning(sprintf("Reference mismatch at position %d: expected %s, got %s",
                      pos, ref_base, ref_in_seq))
    }
    
    # Decide substitution base
    if (allele_idx == 0) {
      base_to_use <- ref_base
    } else if (allele_idx > length(alt_bases)) {
      warning(sprintf("Allele index %d exceeds ALT count at %d", allele_idx, pos))
      next
    } else {
      base_to_use <- alt_bases[allele_idx]
    }
    
    # Replace base(s) in sequence
    hap_seq <- paste0(
      substring(hap_seq, 1, idx - 1),
      base_to_use,
      substring(hap_seq, idx + ref_len)
    )
  }
  
  # ==== [6] Output FASTA File ====
  writeLines(
    c(paste0(">", hap), hap_seq),
    con = paste0("fasta_output/", hap, ".fasta")
  )
}
