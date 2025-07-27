# ==============================================================================
# 7.Fig5BC_Olink_plot.R
# Author: Xiaopu Zhou
# Purpose:
#   - Generates volcano plot for Olink-association results (score_leading only)
#   - Performs GO/KEGG enrichment for significant hits (FDR < 0.05)
#   - Outputs static and interactive enrichment plots
# ==============================================================================

# Set working directory
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure5_UKB/04.Olink/plot")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

data <- read.table("../Olink_association_results_36021.raw.txt", header = TRUE, stringsAsFactors = FALSE)

# Replace extreme FDR for log scale
data$FDR[data$FDR < 1e-300] <- 1e-300

# Focus only on score_leading associations
data <- data[data$Predictor == "score_leading", ]

# ------------------------------------------------------------------------------
# Volcano Plot
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

# Prepare data
volcano_df <- data %>%
  mutate(
    negLog10FDR = -log10(FDR),
    sig = case_when(
      FDR < 0.05 & Beta > 0 ~ "Up",
      FDR < 0.05 & Beta < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Create volcano plot
volcano_plot <- ggplot(volcano_df, aes(x = Beta, y = negLog10FDR)) +
  geom_point(aes(colour = sig), alpha = 1, size = 2) +
  scale_colour_manual(values = c(Up = "#a357e6", Down = "#6ff7aa", NS = "grey80")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = "Beta",
    y = expression(-log[10](FDR)),
    colour = "Significance"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

# Save plot
png("Olink_volcano_leading.png", res = 1200, width = 12, height = 5, pointsize = 10, units = "cm")
volcano_plot
dev.off()

# ------------------------------------------------------------------------------
# GO/KEGG Enrichment for Significant Hits
# ------------------------------------------------------------------------------

library(gprofiler2)

# Get significant protein list
data_fdr <- data[data$FDR < 0.05, ]
gene_symbols <- data_fdr$Region  # assuming protein symbols are stored in Region column

# Run enrichment
ego <- gost(
  query = gene_symbols,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  correction_method = "fdr",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
  user_threshold = 0.05
)

# Convert results to data.frame
ego_df <- if (is.list(ego$result) && !is.data.frame(ego$result)) {
  list_rbind(ego$result, names_to = "query")
} else {
  ego$result
}

# Flatten list columns
ego_df <- ego_df %>%
  mutate(across(where(is.list), ~ sapply(.x, paste, collapse = "; ")))

# Save result table
write.csv(ego_df, "GO_enrichment_Olink_leading_FDR_0.05.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Static Dotplot
# ------------------------------------------------------------------------------

library(ggplot2)

# Non-interactive dotplot
p <- gostplot(ego, interactive = FALSE, capped = FALSE)

# Adjust bubble size and remove legend if desired
p_small <- p +
  scale_size_area(max_size = 1.5, guide = "none")

# Save plot
png("GO_dotplot.png", res = 1200, width = 12, height = 5, pointsize = 10, units = "cm")
p_small
dev.off()

# Optional ggplot2-based ggsave as backup
ggsave("GO_dotplot.png", plot = p_small,
       width = 12, height = 5, units = "cm", dpi = 1200, pointsize = 10)

# ------------------------------------------------------------------------------
# Interactive HTML Dotplot (Optional)
# ------------------------------------------------------------------------------

p_interactive <- gostplot(ego, interactive = TRUE, capped = TRUE)

publish_gostplot(
  p_interactive,
  highlight_terms = ego$result$term_id[1:10],
  filename = "GO_enrichment_interactive.html"
)

# ==============================================================================
# End of script
# ==============================================================================
