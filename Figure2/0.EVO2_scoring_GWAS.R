# ------------------------------------------------------------------------------
# 0.EVO2_scoring_GWAS.R
# Author: Xiaopu Zhou
# Description:
#   This script performs reference allele matching for APOE region SNPs
#   and generates mutant FASTA sequences based on allele mismatches.
#   Outputs include:
#     - Filtered SNPs with valid reference/alt alleles
#     - Ref/alt matched SNPs
#     - Mutated FASTA sequences
# ------------------------------------------------------------------------------

# === Configuration ===

# Coordinates for APOE region (GRCh38)
region_start <- 45345000
region_end   <- 45425000
region_length <- region_end - region_start + 1

# Input files
snp_file <- "SNP_demo.txt"
ref_file <- "APOE_ref_demo.txt"

# Output directory
output_dir <- "results"
fasta_dir  <- file.path(output_dir, "FASTA_modified")

# Create output directories if missing
dir.create(output_dir, showWarnings = FALSE)
dir.create(fasta_dir, showWarnings = FALSE)

# === Load Data ===

message("Loading reference sequence...")
ref_seq <- readLines(ref_file)[1]
stopifnot(nchar(ref_seq) == region_length)

message("Loading SNP data...")
snp <- read.table(snp_file, stringsAsFactors = FALSE, header = TRUE)

# === Helper Functions ===

# Get base from reference sequence
get_ref_base <- function(bp) {
  substr(ref_seq, bp - region_start + 1, bp - region_start + 1)
}

# Replace base at position in sequence
replace_base <- function(seq, pos, new_base) {
  substr(seq, pos, pos) <- new_base
  return(seq)
}

# Write mutated FASTA to file
write_mutated_fasta <- function(snp_row, from_col, to_col, out_dir) {
  bp <- as.numeric(snp_row$BP)
  from_base <- snp_row[[from_col]]
  to_base <- snp_row[[to_col]]
  snp_id <- snp_row$SNP

  pos <- bp - region_start + 1
  mutated_seq <- replace_base(ref_seq, pos, to_base)

  file_name <- paste0(snp_id, "_", from_base, "_to_", to_base, ".txt")
  writeLines(mutated_seq, con = file.path(out_dir, file_name))
}

# === Annotation ===

message("Annotating SNPs against reference sequence...")

# Check if ref_seq base matches A2 or A1 allele
snp$match_ref <- mapply(function(bp, a2) get_ref_base(bp) == a2, snp$BP, snp$A2)
snp$match_alt <- mapply(function(bp, a1) get_ref_base(bp) == a1, snp$BP, snp$A1)

# Keep SNPs with at least one match
snp_clean <- snp[snp$match_ref | snp$match_alt, ]

# Save full matched list
write.table(snp_clean,
            file = file.path(output_dir, "APOE_ref_with_match_demo.txt"),
            row.names = FALSE, quote = FALSE, sep = '\t')

# Split into REF and ALT matches
snp_ref <- snp_clean[snp_clean$match_ref, ]
snp_alt <- snp_clean[snp_clean$match_alt, ]

# Save individual match files
write.table(snp_ref,
            file = file.path(output_dir, "APOE_ref_REF_match_demo.txt"),
            row.names = FALSE, quote = FALSE, sep = '\t')

write.table(snp_alt,
            file = file.path(output_dir, "APOE_ref_ALT_match_demo.txt"),
            row.names = FALSE, quote = FALSE, sep = '\t')

# === Generate Mutated FASTA Sequences ===

message("Generating mutated FASTA sequences...")

# REF match: modify from A2 → A1
apply(snp_ref, 1, function(row) {
  write_mutated_fasta(as.list(row), from_col = "A2", to_col = "A1", out_dir = fasta_dir)
})

# ALT match: modify from A1 → A2
apply(snp_alt, 1, function(row) {
  write_mutated_fasta(as.list(row), from_col = "A1", to_col = "A2", out_dir = fasta_dir)
})

message("All outputs written to '", output_dir, "' and FASTAs to '", fasta_dir, "'")
