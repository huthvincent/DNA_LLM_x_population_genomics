# ==============================================================================
# 0.EVO2_scoring_ADNI.R
# Author: Xiaopu Zhou
# Purpose:
#   - Load cleaned APOE reference sequence from FASTA-like format
#   - Parse individual VCF files from various ADNI datasets
#   - For each sample, reconstruct two phased haplotypes based on genotype
#   - Output per-sample FASTA files for downstream alignment or analysis
# ==============================================================================

# Set working directory

# Load reference sequence
ref <- read.table("../APOE_ref_45345000_45425000_clean.txt", stringsAsFactors = FALSE, header = FALSE)
ref_seq <- ref$V1
stopifnot(nchar(ref_seq) == 45425000 - 45345000 + 1)

# Define dataset names (matching VCF file names in ./APOE_extract/)
datasets <- c(
  "ADNI3_PLINK_Final",
  "ADNI_GO2_GWAS_2nd_orig_BIN",
  "ADNI_GO_2_Forward_Bin",
  "ADNI_cluster_01_forward_757LONI",
  "WGS_Omni25_BIN_wo_ConsentsIssues",
  "ADNI_WGS_808"
)

# Loop over each dataset
for (filename in datasets) {
  
  # Compose file path
  vcf_path <- paste0("./APOE_extract/", filename, ".APOE_211_site.vcf")
  
  # Read VCF (skipping header lines starting with "##")
  skip_lines <- length(grep("^##", readLines(vcf_path)))
  snp_table <- read.table(
    file = vcf_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    skip = skip_lines
  )
  
  # Extract sample names
  sample_names <- colnames(snp_table)[10:ncol(snp_table)]
  
  # Create output directory
  dir.create(paste0("fasta_output/", filename), showWarnings = FALSE, recursive = TRUE)
  
  # Process each sample
  for (sample in sample_names) {
    
    # Duplicate reference for two haplotypes
    hap1 <- ref_seq
    hap2 <- ref_seq
    
    for (i in 1:nrow(snp_table)) {
      pos <- snp_table$POS[i]
      ref_base <- snp_table$REF[i]
      alt_base <- snp_table$ALT[i]
      
      genotype <- snp_table[i, sample]
      genotype <- strsplit(genotype, ":")[[1]][1]
      alleles <- strsplit(genotype, "\\|")[[1]]
      
      idx <- pos - 45345000 + 1  # convert genomic position to 1-based index
      
      # Update hap1 and hap2 based on genotype
      if (alleles[1] == "1") {
        substr(hap1, idx, idx) <- alt_base
      } else {
        substr(hap1, idx, idx) <- ref_base
      }
      
      if (alleles[2] == "1") {
        substr(hap2, idx, idx) <- alt_base
      } else {
        substr(hap2, idx, idx) <- ref_base
      }
    }
    
    # Write to FASTA
    writeLines(
      c(paste0(">", sample, "_hap1"), hap1),
      con = paste0("fasta_output/", filename, "/", filename, "-", sample, "_hap1.fasta")
    )
    writeLines(
      c(paste0(">", sample, "_hap2"), hap2),
      con = paste0("fasta_output/", filename, "/", filename, "-", sample, "_hap2.fasta")
    )
  }
}
