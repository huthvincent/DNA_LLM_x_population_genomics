# ------------------------------------------------------------------------------
# 2.Fig2B.R
# Author: Xiaopu Zhou
# Description:
#   Generate two separate scatter plots:
#     (1) EVO2 score (Delta) vs. Z-score
#     (2) EVO2 score (Delta) vs. GWAS effect size (BETA)
#   Each with linear fit, 95% CI, and axis guides.
# ------------------------------------------------------------------------------


# === Load Required Libraries ===
library(ggplot2)
library(ggsci)
library(dplyr)

# === Load Input Data ===
gwas_data <- read.table("GWAS_delta_summary_211.txt", header = TRUE, stringsAsFactors = FALSE)

# Optional: Load clumped SNPs list (not used directly here)
clumped_set <- read.table("Clumped_28_list.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# === Set Colors ===
point_color <- "#3B4992FF"   # AAAS blue
smooth_fill <- "#A6BCDFFF"   # lighter AAAS blue for confidence interval

# ------------------------------------------------------------------------------
# Plot 1: EVO2 Delta vs. Z-score
# ------------------------------------------------------------------------------

png(file = "EVO2_Z_value_AllSite.png", res = 1200, width = 9, height = 9, pointsize = 10, units = "cm")

ggplot(gwas_data, aes(x = Delta, y = Z)) +
  geom_point(color = point_color, size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              color = point_color, fill = smooth_fill,
              size = 1.2, alpha = 0.4) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "EVO2 Score vs. Z-score",
    x = "EVO2 Score (Delta)",
    y = "GWAS Z-score"
  ) +
  xlim(-30, 20) +
  ylim(-40, 60) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2)

dev.off()

# ------------------------------------------------------------------------------
# Plot 2: EVO2 Delta vs. GWAS Effect Size (BETA)
# ------------------------------------------------------------------------------

png(file = "EVO2_effect_size_AllSite.png", res = 1200, width = 9, height = 9, pointsize = 10, units = "cm")

ggplot(gwas_data, aes(x = Delta, y = BETA)) +
  geom_point(color = point_color, size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              color = point_color, fill = smooth_fill,
              size = 1.2, alpha = 0.4) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "EVO2 Score vs. Effect Size",
    x = "EVO2 Score (Delta)",
    y = "GWAS Effect Size (BETA)"
  ) +
  xlim(-30, 20) +
  ylim(-0.2, 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2)

dev.off()
