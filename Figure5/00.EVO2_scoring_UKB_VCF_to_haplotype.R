# ==============================================================================
# 00.EVO2_scoring_UKB_VCF_to_haplotype.R
# Author: Xiaopu Zhou
# Purpose:
#   - Read WGS VCF containing phased APOE-region variants (UK Biobank subset)
#   - Extract genotype calls per haplotype (0/1 format)
#   - Output long-format haplotype table per individual and allele copy
# ==============================================================================

# ==== [0] Set Working Directory ====

# ==== [1] Load Required Libraries ====
library(vcfR)
library(data.table)

# ==== [2] Load VCF File ====
vcf_path <- "./VCF/UKB_WGS_APOE.AD_GWAS_211_site.Final.vcf"
vcf <- read.vcfR(vcf_path)

# ==== [3] Extract Variant Metadata ====
fix_df <- as.data.frame(getFIX(vcf))
fix_df$POS <- as.numeric(fix_df$POS)

# Optional: filter by genomic range (GRCh38 APOE locus window)
pos_min <- 44841743
pos_max <- 44921743
vcf_filtered <- vcf[fix_df$POS >= pos_min & fix_df$POS <= pos_max, ]
fix_df <- as.data.frame(getFIX(vcf_filtered))  # re-extract after filtering
fix_df$POS <- as.numeric(fix_df$POS)

# ==== [4] Extract Genotypes ====
gt <- extract.gt(vcf_filtered, element = "GT")

# Convert to data.frame and add rownames as variant ID
gt_df <- as.data.frame(gt)
gt_df$VARIANT <- rownames(gt_df)

# ==== [5] Split GT to Hap1 and Hap2 ====
hap_split <- function(x, which_hap = 1) {
  sapply(strsplit(x, "\\|"), function(y) as.integer(y[which_hap]))
}

# Drop "VARIANT" column for haplotype processing
hap1_df <- as.data.frame(lapply(gt_df[, -ncol(gt_df)], hap_split, which_hap = 1))
hap2_df <- as.data.frame(lapply(gt_df[, -ncol(gt_df)], hap_split, which_hap = 2))

# Rename Columns
colnames(hap1_df) <- paste0(colnames(hap1_df), "_Hap1")
colnames(hap2_df) <- paste0(colnames(hap2_df), "_Hap2")

# ==== [6] Combine: Variant Info + Haplotypes ====
variant_info <- fix_df[, c("CHROM", "POS", "ID", "REF", "ALT")]
final_df <- cbind(variant_info, hap1_df, hap2_df)

# ==== [7] Save Output ====
write.table(
  final_df,
  file = "./All_site_genotype/UKB_490558_203_haplotype_table.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==== [8] Preview ====
head(final_df)
