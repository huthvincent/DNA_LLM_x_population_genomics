# ==============================================================================
# 5.Fig3E.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merges pangenome delta scores with sample metadata
#   - Filters for APOE4 haplotypes
#   - Visualizes delta distribution across ancestries using ggridges
#   - Outputs summary statistics and saves plot
# ==============================================================================

# ==== [1] Setup ====

# ==== [2] Load required libraries ====
library(ggplot2)
library(ggridges)
library(dplyr)

# ==== [3] Load input data ====
# EVO2 delta score per haplotype
data <- read.table("Pangenome_EVO2_score.txt", header = TRUE, stringsAsFactors = FALSE)

# Standardize APOE SNP column names
names(data)[names(data) %in% c("rs429358", "rs7412")] <- c("APOE4", "APOE2")

# Sample-level metadata
anno <- read.table("igsr_samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge data with annotation
data <- merge(data, anno, by.x = "IID", by.y = "Sample_name")

# ==== [4] Derive APOE genotype label ====
data$APOE3 <- 1 - data$APOE4 - data$APOE2

data$APOE <- NA
data$APOE[data$APOE4 == 1] <- "APOE4"
data$APOE[data$APOE2 == 1] <- "APOE2"
data$APOE[data$APOE3 == 1] <- "APOE3"

data$APOE <- relevel(as.factor(data$APOE), ref = "APOE3")

# ==== [5] Subset to APOE4 haplotypes ====
data_APOE4 <- data[data$APOE4 == 1, ]

# Set ancestry display order
data_APOE4$Superpopulation_code <- factor(
  data_APOE4$Superpopulation_code,
  levels = rev(c("AFR", "AMR", "EAS", "EUR", "SAS"))
)

# ==== [6] Define colors ====
npg_color_lines <- c(
  AFR = "#B22222",  # Firebrick
  AMR = "#1E88A8",  # Blue
  EAS = "#00775E",  # Teal
  EUR = "#2B3E68",  # Indigo
  SAS = "#C46245"   # Burnt sienna
)

# ==== [7] Plot ridge plot ====
plot <- ggplot(data_APOE4, aes(x = delta, y = Superpopulation_code, fill = Superpopulation_code)) +
  geom_density_ridges(scale = 0.8, alpha = 0.8) +
  scale_fill_manual(values = npg_color_lines) +
  theme_classic() +
  xlim(-0.002, 0.001) +
  theme(legend.position = "none")

# Save figure
png("APOE4_all_ancestry.png", res = 1200, width = 12, height = 6, pointsize = 10, units = "cm")
print(plot)
dev.off()

# Output full merged table for reproducibility
write.table(data, "Pangenome_APOE4_hap_score_plot_Final.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ==== [8] Summary table per ancestry ====
summary_data <- data_APOE4 %>%
  group_by(Superpopulation_code) %>%
  summarise(
    mean_EVO2 = mean(delta, na.rm = TRUE),
    sem_EVO2 = sd(delta, na.rm = TRUE) / sqrt(sum(!is.na(delta)))
  )

print(summary_data)
