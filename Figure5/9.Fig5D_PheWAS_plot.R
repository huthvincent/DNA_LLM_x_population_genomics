# ==============================================================================
# 9.Fig5D_PheWAS_plot.R
# Author: Xiaopu Zhou
# Purpose:
#   - Generate Manhattan-style plot for PheWAS results from EVO2 score_leading
#   - Color-coded by Phecode category, highlighting FDR-significant traits
#   - Output: PNG figure of plot
# ==============================================================================

# ------------------------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------------------------
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure5_UKB/05.PheWAS")

# ------------------------------------------------------------------------------
# Load association results and phecode annotation
# ------------------------------------------------------------------------------
data <- read.table("PheWAS_association_results_1547disease_n331925.raw.txt",
                   header = TRUE, stringsAsFactors = FALSE)

anno <- read.csv("Phecode_map12_filtered.csv", header = TRUE, stringsAsFactors = FALSE)
anno <- anno[, c("Phecode", "PhecodeString", "PhecodeCategory")]
anno <- anno[!duplicated(anno), ]

data <- merge(data, anno, by.x = "Disease", by.y = "Phecode")
data <- data[!is.na(data$PhecodeCategory), ]

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)       # AAAS color palette
library(colorspace)  # for lighten()

# ------------------------------------------------------------------------------
# Filter for score_leading and calculate FDR
# ------------------------------------------------------------------------------
phewas_df <- data[data$Predictor == "score_leading", ]
phewas_df$FDR <- p.adjust(phewas_df$p_value, method = "fdr")

# ------------------------------------------------------------------------------
# Prepare plotting data
# ------------------------------------------------------------------------------
plot_df <- phewas_df %>%
  mutate(
    logP = -log10(p_value),
    PhecodeCategory = factor(PhecodeCategory, levels = unique(PhecodeCategory)),
    sig = FDR < 0.05
  ) %>%
  arrange(PhecodeCategory, p_value) %>%
  group_by(PhecodeCategory) %>%
  mutate(order_in_cat = row_number()) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Compute X-axis category offsets
# ------------------------------------------------------------------------------
cat_sizes <- plot_df %>%
  count(PhecodeCategory, name = "n_in_cat") %>%
  mutate(offset = lag(cumsum(n_in_cat), default = 0))

plot_df <- plot_df %>%
  left_join(cat_sizes, by = "PhecodeCategory") %>%
  mutate(x_pos = order_in_cat + offset)

tick_df <- cat_sizes %>%
  mutate(tick = offset + n_in_cat / 2)

# ------------------------------------------------------------------------------
# Significance guide and top hits for optional labeling
# ------------------------------------------------------------------------------
sig_cut <- 3e-4
line_y <- -log10(sig_cut)

top_df <- plot_df %>%
  filter(sig) %>%
  group_by(PhecodeCategory) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Construct extended color palette
# ------------------------------------------------------------------------------
n_cat <- nlevels(plot_df$PhecodeCategory)
aaas_dark <- pal_aaas("default")(10)
aaas_light <- lighten(aaas_dark, amount = 0.35)
ext_pal <- as.vector(rbind(aaas_dark, aaas_light))

while (length(ext_pal) < n_cat) {
  ext_pal <- c(ext_pal, as.vector(rbind(lighten(aaas_dark, 0.6), aaas_dark)))
}

reorder_idx <- c(rbind(seq_len(n_cat), rev(seq_len(n_cat))))[1:n_cat]
cat_cols <- ext_pal[reorder_idx]

# ------------------------------------------------------------------------------
# Build plot
# ------------------------------------------------------------------------------
leading_plot <- ggplot(plot_df,
                       aes(x = x_pos, y = logP, colour = PhecodeCategory, size = sig, alpha = sig)) +
  geom_point() +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1.0), guide = "none") +
  scale_colour_manual(values = cat_cols, name = "Phecode category") +
  geom_hline(yintercept = line_y, linetype = "dashed", colour = "grey60") +
  # Optional: Add labels for top phenotype per category
  # geom_text_repel(data = top_df, aes(label = PhecodeString),
  #                 size = 3, segment.color = "grey60", seed = 1) +
  scale_y_continuous(name = expression(-log[10](italic(p))),
                     expand = c(0.01, 0), limits = c(0, 60)) +
  scale_x_continuous(name = "Phecode category",
                     breaks = tick_df$tick,
                     labels = tick_df$PhecodeCategory,
                     expand = c(0.01, 0)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# Save plot
# ------------------------------------------------------------------------------
png(file = "leading_PhwWAS_final.png", res = 1200,
    width = 28, height = 9, pointsize = 10, units = "cm")
plot(leading_plot)
dev.off()

# ==============================================================================
# End of script
# ==============================================================================
