# ==============================================================================
# 4.Fig3D.R
# Author: Xiaopu Zhou
# Purpose:
#   - Classify haplotypes by APOE allele based on rs429358 and rs7412
#   - Compute APOE3 carrier status
#   - Visualize Evo2 delta score distribution per APOE class
# ==============================================================================

# ==== [0] Set Working Directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure3_Pangenome/panel_d")

# ==== [1] Load Required Libraries ====
library(ggplot2)
library(ggridges)

# ==== [2] Load Evo2 Score File with rs429358 & rs7412 ====
data <- read.table("Pangenome_EVO2_score.txt", header = TRUE, stringsAsFactors = FALSE)

# Rename known columns to standard APOE allele labels
names(data)[names(data) %in% c("rs429358", "rs7412")] <- c("APOE4", "APOE2")

# ==== [3] Compute APOE3 (Residual) Haplotypes ====
# APOE3 defined as: neither APOE4 nor APOE2
data$APOE3 <- 1 - data$APOE4 - data$APOE2

# Check allele distribution
table(data$APOE4, data$APOE2, data$APOE3)

# ==== [4] Assign APOE Genotype Label Per Haplotype ====
data$APOE <- NA
data$APOE[data$APOE4 == 1] <- "APOE4"
data$APOE[data$APOE2 == 1] <- "APOE2"
data$APOE[data$APOE3 == 1] <- "APOE3"

# Relevel factor so APOE3 is baseline
data$APOE <- factor(data$APOE, levels = c("APOE3", "APOE2", "APOE4"))

# ==== [5] Ridge Plot: Delta Score by APOE Type ====

ridge_plot <- ggplot(data, aes(x = delta, y = APOE, fill = APOE)) +
  geom_density_ridges(scale = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("grey80", "#7bdb90", "#a85bc9")) +
  theme_classic(base_size = 12) +
  xlim(-0.002, 0.001) +
  labs(x = "Delta (Evo2 score)", y = "APOE allele") +
  theme(legend.position = "none")

# ==== [6] Save Plot ====
png("APOE234_all_ancestry.png", res = 1200, width = 12, height = 6, pointsize = 10, units = "cm")
print(ridge_plot)
dev.off()

# ==== [7] Output Clean Table for Record ====
write.table(
  data,
  file = "Pangenome_APOE234_hap_score_plot_Final.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
