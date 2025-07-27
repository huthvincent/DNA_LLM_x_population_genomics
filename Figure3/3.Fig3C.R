# ==============================================================================
# 3.Fig3C.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merge RFMix2 ancestry output with Evo2 delta scores
#   - Pivot region-ancestry data
#   - Compute partial Spearman correlation for each Region × Ancestry
#   - Output correlation and p-value matrices
# ==============================================================================


# ==== [1] Load Required Libraries ====
library(dplyr)
library(tidyr)
library(ppcor)

# ==== [2] Read and Merge Data ====

# RFMix2 output (per haplotype, per region ancestry)
rfmix_data <- read.table("Pangenome.RFMIX2.clean.txt", header = TRUE, stringsAsFactors = FALSE)

# Evo2 delta score per haplotype
evo2_data <- read.table("Pangenome_EVO2.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge datasets on haplotype ID
merged_data <- merge(rfmix_data, evo2_data, by.x = "ID", by.y = "Hap_ID")

# Save merged data (optional)
write.table(merged_data, "Pangenome.RFMIX2.clean.EVO2.txt", quote = FALSE, sep = '\t', row.names = FALSE)

# ==== [3] Reshape to Region × Ancestry Wide Format ====

long_data <- merged_data %>%
  pivot_longer(
    cols = starts_with("X19"),  # local ancestry region columns
    names_to = "Region",
    values_to = "Value"
  ) %>%
  mutate(Region_Ancestry = paste(Region, Ethnicity, sep = "_")) %>%
  select(ID, delta, Region_Ancestry, Value) %>%
  pivot_wider(names_from = Region_Ancestry, values_from = Value)

# ==== [4] Partial Correlation by Region × Ancestry ====

# Identify regions and ancestry groups
region_cols <- grep("^X19", names(merged_data), value = TRUE)
ancestries <- unique(merged_data$Ethnicity)

# Initialize results dataframe
results <- data.frame()

for (region in region_cols) {
  ancestry_cols <- paste0(region, "_", ancestries)
  ancestry_cols <- ancestry_cols[ancestry_cols %in% names(long_data)]  # keep valid cols
  
  for (anc in ancestries) {
    target <- paste0(region, "_", anc)
    controls <- setdiff(ancestry_cols, target)
    
    if (target %in% names(long_data) && length(controls) > 0) {
      df <- long_data[, c("delta", target, controls)]
      df <- na.omit(df)
      
      # At least 6 samples and more than 1 unique value to compute correlation
      if (nrow(df) > 5 && length(unique(df[[target]])) > 1) {
        pc <- pcor.test(x = df[[target]], y = df$delta, z = df[, controls[-1]], method = "spearman")
        
        results <- rbind(results, data.frame(
          Region = region,
          Ancestry = anc,
          Correlation = pc$estimate,
          P_value = pc$p.value
        ))
      }
    }
  }
}

# ==== [5] Reshape Results to Correlation and P-value Matrices ====

cor_matrix <- results %>%
  select(Region, Ancestry, Correlation) %>%
  pivot_wider(names_from = Ancestry, values_from = Correlation)

p_matrix <- results %>%
  select(Region, Ancestry, P_value) %>%
  pivot_wider(names_from = Ancestry, values_from = P_value)

# ==== [6] Output Results ====

# Print to console
print("=== Partial Correlation Matrix ===")
print(cor_matrix)

print("=== P-value Matrix ===")
print(p_matrix)

# Optionally save results
write.table(cor_matrix, "RegionAncestry_PartialCorrelation.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(p_matrix,   "RegionAncestry_Pvalues.txt",         sep = "\t", quote = FALSE, row.names = FALSE)
