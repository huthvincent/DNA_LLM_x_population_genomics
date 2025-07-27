# ==============================================================================
# 5.Fig5A_SWI_plot.R
# Author: Xiaopu Zhou
# Purpose:
#   - Visualizes T2star and susceptibility-weighted imaging (SWI) metrics
#     using cerebroViz
#   - Separate plots are created for score_leading and score_alternative
#     on both left and right hemispheres
#   - Results are saved into distinct output directories for each data type
# ==============================================================================

# Load libraries
library(cerebroViz)

# ========== Block 1: T2star Plotting ==========

# Set working directory for T2star data

data_t2 <- read.table("T2star_plot.txt", header = TRUE, stringsAsFactors = FALSE)

# Create output directory
dir.create("T2star_plot/", showWarnings = FALSE)

# Function to generate cerebroViz plot
plot_cerebro <- function(data, side, predictor, output_prefix) {
  dt <- subset(data, side == !!side & Predictor == !!predictor)
  dt <- dt[order(dt$p_value), ]
  dt <- dt[!duplicated(dt$Region_cerebroViz), c("Region_cerebroViz", "t_value")]
  
  mat <- as.matrix(dt[, "t_value"])
  rownames(mat) <- dt$Region_cerebroViz
  colnames(mat) <- predictor
  attributes(mat)$class <- "matrix"
  
  # Add dummy values to fix color scale
  mat_range <- rbind(mat, .min_dummy = -4, .max_dummy = 4)
  attributes(mat_range)$class <- "matrix"
  
  cerebroViz(mat_range,
             filePrefix = output_prefix,
             timePoint = 1,
             naHatch = FALSE,
             regLabel = TRUE,
             palette = c("#6ff7aa", "white", "#a357e6"),
             secPalette = c("grey80", "black", "white"))
}

# T2star plots
plot_cerebro(data_t2, "right", "score_leading", "T2star_plot/T2star_right_leading_")
plot_cerebro(data_t2, "left",  "score_leading", "T2star_plot/T2star_left_leading_")
plot_cerebro(data_t2, "left",  "score_alternative", "T2star_plot/T2star_left_alternative_")
plot_cerebro(data_t2, "right", "score_alternative", "T2star_plot/T2star_right_alternative_")


# ========== Block 2: Susceptibility Plotting ==========

# Load susceptibility data
data_susc <- read.table("susceptibility_plot.txt", header = TRUE, stringsAsFactors = FALSE)

# Create output directory
dir.create("susceptibility_plot/", showWarnings = FALSE)

# Susceptibility plots
plot_cerebro(data_susc, "right", "score_leading", "susceptibility_plot/susceptibility_right_leading_")
plot_cerebro(data_susc, "left",  "score_leading", "susceptibility_plot/susceptibility_left_leading_")
plot_cerebro(data_susc, "left",  "score_alternative", "susceptibility_plot/susceptibility_left_alternative_")
plot_cerebro(data_susc, "right", "score_alternative", "susceptibility_plot/susceptibility_right_alternative_")

# ==============================================================================
# End of script
# ==============================================================================
