# ------------------------------------------------------------------------------
# Evo2 LDSC PheWAS plot
# Author: Xiaopu Zhou
# Description:
# This script plots LDSC coefficient z-scores across UKB phecode categories.
# ------------------------------------------------------------------------------


# === Load Required Libraries ===

library(readr)
library(dplyr)
library(ggplot2)

# === Load Data ===

data <- read_tsv(
  "Evo2_LDSC_UKB.txt",
  show_col_types = FALSE
)

# === Define Colors ===

category_colors <- c(
  "circulatory system" = "#0072B2",
  "endocrine/metabolic" = "#009E73",
  "mental disorders" = "#56B4E9",
  "digestive" = "#E69F00",
  "musculoskeletal" = "#D55E00",
  "respiratory" = "#00BFC4",
  "neurological" = "#332288",
  "sense organs" = "#CC79A7",
  "genitourinary" = "#44AA99",
  "dermatologic" = "#999933",
  "infectious diseases" = "#F0E442",
  "hematopoietic" = "#AA4499",
  "pregnancy complications" = "#882255",
  "congenital anomalies" = "#BBBBBB",
  "injuries & poisonings" = "#666666",
  "neoplasms" = "#CC3311"
)

# === Prepare Data ===

phewas_df <- data %>%
  mutate(
    bonferroni = p.adjust(Enrichment_p, method = "bonferroni"),
    z = Coefficient_z_score
  )

z_cut <- phewas_df %>%
  filter(bonferroni < 0.05, is.finite(z)) %>%
  summarise(z_cut = min(abs(z))) %>%
  pull(z_cut)

phewas_df <- phewas_df %>%
  mutate(sig = abs(z) >= z_cut)

write.table(
  phewas_df,
  "Evo2_LDSC_UKB_bonferroni.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# === Order Categories ===

cat_levels <- phewas_df %>%
  group_by(phecode_category_Final) %>%
  summarise(sentinel_z = max(z, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(sentinel_z)) %>%
  pull(phecode_category_Final)

plot_df <- phewas_df %>%
  mutate(
    phecode_category_Final = trimws(as.character(phecode_category_Final)),
    phecode_category_Final = factor(phecode_category_Final, levels = cat_levels)
  ) %>%
  arrange(phecode_category_Final, desc(z)) %>%
  group_by(phecode_category_Final) %>%
  mutate(order_in_cat = row_number()) %>%
  ungroup()

# === Compute X-axis Positions ===

cat_sizes <- plot_df %>%
  count(phecode_category_Final, name = "n_in_cat") %>%
  mutate(offset = lag(cumsum(n_in_cat), default = 0))

plot_df <- plot_df %>%
  left_join(cat_sizes, by = "phecode_category_Final") %>%
  mutate(x_pos = order_in_cat + offset)

tick_df <- cat_sizes %>%
  mutate(tick = offset + n_in_cat / 2)

# === Generate Plot ===

leading_plot <- ggplot(
  plot_df,
  aes(
    x = x_pos,
    y = z,
    colour = phecode_category_Final,
    size = sig,
    alpha = sig
  )
) +
  geom_point() +
  geom_hline(
    yintercept = z_cut,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_colour_manual(
    values = category_colors,
    name = "Phecode category"
  ) +
  scale_size_manual(
    values = c(`FALSE` = 1, `TRUE` = 2),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(`FALSE` = 0.3, `TRUE` = 1),
    guide = "none"
  ) +
  scale_y_continuous(
    name = "Coefficient z-score",
    limits = c(-5, 24),
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(
    name = "Phecode category",
    breaks = tick_df$tick,
    labels = tick_df$phecode_category_Final,
    expand = c(0.01, 0)
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# === Save Plot ===

png(
  file = "leading_PheWAS_zscore_ranked_empirical_bonf.png",
  res = 1200,
  width = 19,
  height = 9,
  pointsize = 5,
  units = "cm"
)

plot(leading_plot)

dev.off()

message("Empirical Bonferroni z cutoff = ", round(z_cut, 3))