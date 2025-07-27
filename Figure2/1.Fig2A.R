# ------------------------------------------------------------------------------
# 1.Fig2A.R
# Author: Xiaopu Zhou
# Description:
#   This script generates two separate LocusZoom-style plots using `topr`:
#   (1) AD GWAS (lead SNP: rs429358)
#   (2) EVO2 score (lead SNP: rs7412)
# ------------------------------------------------------------------------------

# === Set Working Directory ===
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure2_GWAS/panel_a")

# === Load Required Libraries ===
if (!requireNamespace("topr", quietly = TRUE)) {
  devtools::install_github("totajuliusd/topr")
}
library(topr)
library(dplyr)
library(LDlinkR)
library(stringr)

# === Load GWAS Summary Data ===
gwas_data <- read.table("GWAS_delta_summary_211.txt", header = TRUE, stringsAsFactors = FALSE)

# === Helper Function: Prepare LD & Merge with GWAS ===
prepare_ld_data <- function(df, lead_snp_id, p_col_index, label) {
  # Rename columns for compatibility
  df <- df[, c(9, 10, 2, 3, 4, p_col_index, 21, 22)]
  colnames(df) <- c("CHROM", "POS", "ID", "A1", "A2", "P", "SYMBOL", "Consequence")
  
  # Get SNP coordinates
  lead_snp <- df %>% filter(ID == lead_snp_id)
  snp_coord <- paste0("chr", lead_snp$CHROM, ":", lead_snp$POS)
  
  # Query LD data using LDlinkR
  LD_proxy_data <- LDproxy(
    snp = snp_coord,
    pop = "EUR",
    r2d = "r2",
    token = "1ea30fd99d17",  # Replace with your valid token
    genome_build = "grch37"
  )
  
  # Format LD and merge
  LD_info <- LD_proxy_data %>%
    mutate(POS = as.numeric(str_extract(Coord, "(?<=:)[0-9]+"))) %>%
    select(POS, R2)
  
  merged <- merge(df, LD_info, by = "POS", all.x = TRUE)
  merged$R2[is.na(merged$R2)] <- 0
  return(merged)
}

# === Prepare Data for Each Plot ===

# 1. AD GWAS (P-value in column 15)
snps_ld_ad <- prepare_ld_data(gwas_data, lead_snp_id = "rs429358", p_col_index = 15, label = "AD GWAS")

# 2. EVO2 Score (P-value in column 6)
snps_ld_evo2 <- prepare_ld_data(gwas_data, lead_snp_id = "rs7412", p_col_index = 6, label = "EVO2 Score")

# === Plot 1: AD GWAS ===
png("AD_GWAS_locuszoom.png", res = 1200, width = 16, height = 8, units = "cm", pointsize = 10)

locuszoom(
  snps_ld_ad,
  region = "chr19:45345000-45425000",
  build = "GRCh37",
  log_trans_p = FALSE,
  rsids = c("rs429358"),
  label_fontface = "italic",
  gene_color = "#113082",
  unit_gene = 1.5,
  protein_coding_only = TRUE,
  show_exons = TRUE,
  sign_thresh = -1,
  gene_label_size = 0,
  ymax = 60,
  ymin = 0,
  size = 1.3
)

dev.off()

# === Plot 2: EVO2 Score ===
png("EVO2_locuszoom.png", res = 1200, width = 16, height = 8, units = "cm", pointsize = 10)

locuszoom(
  snps_ld_evo2,
  region = "chr19:45345000-45425000",
  build = "GRCh37",
  log_trans_p = FALSE,
  rsids = c("rs7412", "rs429358"),
  label_fontface = "italic",
  gene_color = "#113082",
  unit_gene = 1.5,
  protein_coding_only = TRUE,
  show_exons = TRUE,
  sign_thresh = -1,
  gene_label_size = 0,
  ymax = 30,
  ymin = 0,
  size = 1.3
)

dev.off()
