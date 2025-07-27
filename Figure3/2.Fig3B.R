# ==============================================================================
# 2.Fig3B.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merge per-haplotype APOE expression and delta score
#   - Annotate haplotypes as leading vs. alternative
#   - Visualize correlation between delta and expression
#   - Fit separate linear mixed models by haplotype type
# ==============================================================================


# ==== [1] Load Libraries ====
library(dplyr)
library(ggplot2)
library(robustbase)
library(lme4)       # For linear mixed models

# ==== [2] Read Input Data ====

# APOE normalized expression data for 33 samples (per haplotype)
APOE <- read.table("Pangenome_APOE_33_expression.txt", header = TRUE, stringsAsFactors = FALSE)

# Per-haplotype delta score (EVO2)
EVO2 <- read.table("Pangenome_EVO2_score.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge by haplotype/sample ID
APOE_anno <- merge(APOE, EVO2, by.x = "ID_clean", by.y = "IID")

# ==== [3] Annotate Haplotype Rank: Leading vs Alternative ====

APOE_anno <- APOE_anno %>%
  group_by(ID_clean) %>%
  arrange(desc(delta), .by_group = TRUE) %>%
  mutate(rank_label = if_else(row_number() == 1, "leading", "alternative")) %>%
  ungroup()

# ==== [4] Define Colors ====
custom_colors <- c(
  leading = "#984ea3",      # Purple
  alternative = "#66c2a5"   # Light green
)

# ==== [5] Generate Scatter Plot ====

png("APOE_indiviudal_expression_correlation.png", res = 1200,
    width = 13, height = 9, pointsize = 10, units = "cm")

ggplot(APOE_anno, aes(x = delta, y = APOE_norm)) +
  geom_point(aes(color = rank_label), size = 1.5, alpha = 0.9) +
  geom_smooth(
    data = subset(APOE_anno, rank_label == "leading"),
    aes(color = rank_label, fill = rank_label),
    method = "lm", se = TRUE, fullrange = TRUE, size = 1, alpha = 0.15
  ) +
  geom_smooth(
    data = subset(APOE_anno, rank_label == "alternative"),
    aes(color = rank_label, fill = rank_label),
    method = "lm", se = TRUE, fullrange = TRUE, size = 1, alpha = 0.15
  ) +
  scale_color_manual(values = custom_colors, name = "Haplotype Type") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  theme_classic(base_size = 12) +
  labs(x = "Delta (per haplotype)", y = "Normalized APOE Expression")

dev.off()

# ==== [6] Fit Linear Mixed Models by Haplotype Type ====

# Fit model: APOE_norm ~ delta + (1|Sex) for leading haplotypes
leading_model <- lmer(APOE_norm ~ delta + (1 | Sex), data = subset(APOE_anno, rank_label == "leading"))

# Fit model: APOE_norm ~ delta + (1|Sex) for alternative haplotypes
alternative_model <- lmer(APOE_norm ~ delta + (1 | Sex), data = subset(APOE_anno, rank_label == "alternative"))

# ==== [7] Report Model Summaries ====
cat("=== Leading Haplotype Model ===\n")
print(summary(leading_model))
cat("\n=== Confidence Intervals ===\n")
print(confint(leading_model))

cat("\n\n=== Alternative Haplotype Model ===\n")
print(summary(alternative_model))
cat("\n=== Confidence Intervals ===\n")
print(confint(alternative_model))
