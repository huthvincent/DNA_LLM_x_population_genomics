# ------------------------------------------------------------------------------
# 3.Fig2C.R
# Author: Xiaopu Zhou
# Description:
#   Generate two separate correlation plots for 28 clumped SNPs:
#     (1) EVO2 score (Delta) vs. Z-score
#     (2) EVO2 score (Delta) vs. GWAS effect size (BETA)
#   Each includes linear regression fit and 95% CI band.
# ------------------------------------------------------------------------------

# === Set Working Directory ===
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure2_GWAS/Panel_bc")

# === Load Required Libraries ===
library(ggplot2)
library(ggsci)
library(dplyr)

# === Load Input Data ===
gwas_data <- read.table("GWAS_delta_summary_211.txt", header = TRUE, stringsAsFactors = FALSE)
clumped_set <- read.table("Clumped_28_list.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# === Filter for 28 Clumped SNPs ===
gwas_data <- gwas_data %>% filter(rsID %in% clumped_set)

# === Define Plot Parameters ===
main_color <- "#631879"   # Deep purple
fill_color <- "#9755ab"   # Lighter purple for CI
base_font  <- 10

# ------------------------------------------------------------------------------
# Plot 1: EVO2 Delta vs. Z-score
# ------------------------------------------------------------------------------

png("EVO2_Z_value_clumped_28.png", res = 1200, width = 9, height = 9, pointsize = base_font, units = "cm")

ggplot(gwas_data, aes(x = Delta, y = Z)) +
  geom_point(color = main_color, size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              color = main_color, fill = fill_color,
              size = 1.2, alpha = 0.4) +
  theme_classic(base_size = base_font) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(color = "black")
  ) +
  labs(
    title = "EVO2 Score vs. Z-score (Clumped SNPs)",
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

png("EVO2_effect_size_clumped_28.png", res = 1200, width = 9, height = 9, pointsize = base_font, units = "cm")

ggplot(gwas_data, aes(x = Delta, y = BETA)) +
  geom_point(color = main_color, size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              color = main_color, fill = fill_color,
              size = 1.2, alpha = 0.4) +
  theme_classic(base_size = base_font) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(color = "black")
  ) +
  labs(
    title = "EVO2 Score vs. Effect Size (Clumped SNPs)",
    x = "EVO2 Score (Delta)",
    y = "GWAS Effect Size (BETA)"
  ) +
  xlim(-30, 20) +
  ylim(-0.2, 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.2)

dev.off()
