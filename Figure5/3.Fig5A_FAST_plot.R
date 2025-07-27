# ==============================================================================
# 3.Fig5A_FAST_plot.R
# Author: Xiaopu Zhou
# Purpose:
#   - Annotate brain regions for MRI-based FAST association results
#   - Generate cerebroViz plots stratified by hemisphere and score type
#   - Uses t-values from FAST_plot.txt and fixed colormap scale (-4 to 4)
# ==============================================================================

# Set working directory

# Load dependencies
library(dplyr)
library(cerebroViz)

# Load association result
data <- read.table("FAST_plot.txt", stringsAsFactors = FALSE, header = TRUE, sep = '\t')

# Define region-to-macro label mapping (manual curation)
region_map <- c(
  # Limbic
  "Amygdala" = "AMY",
  "Hippocampus" = "HIP",
  "Parahippocampal.Gyrus.anterior.division" = "HIP",
  "Parahippocampal.Gyrus.posterior.division" = "HIP",
  
  # Basal ganglia
  "Caudate" = "CAU",
  "Putamen" = "PUT",
  "Pallidum" = "PUT",
  "Ventral.Striatum" = "STR",
  
  # Thalamus & brainstem
  "Thalamus" = "THA",
  "Brain_Stem" = "BS",
  
  # Cerebellum (collapsed)
  "Crus.I.Cerebellum" = "CB",
  "Crus.II.Cerebellum" = "CB",
  "I.IV.Cerebellum" = "CB",
  "V.Cerebellum" = "CB",
  "VI.Cerebellum" = "CB",
  "VIIb.Cerebellum" = "CB",
  "VIIIa.Cerebellum" = "CB",
  "VIIIb.Cerebellum" = "CB",
  "IX.Cerebellum" = "CB",
  "X.Cerebellum" = "CB",
  
  # Parietal
  "Angular.Gyrus" = "ANG",
  "Supramarginal.Gyrus.anterior.division" = "IPC",
  "Supramarginal.Gyrus.posterior.division" = "IPC",
  "Parietal.Operculum.Cortex" = "IPC",
  "Superior.Parietal.Lobule" = "PL",
  "Precuneous.Cortex" = "PL",
  
  # Sensorimotor
  "Postcentral.Gyrus" = "S1C",
  "Precentral.Gyrus" = "M1C",
  "Juxtapositional.Lobule.Cortex.formerly.Supplementary.Motor.Cortex" = "DFC",
  
  # Occipital / visual
  "Lingual.Gyrus" = "OL",
  "Cuneal.Cortex" = "V1C",
  "Intracalcarine.Cortex" = "V1C",
  "Supracalcarine.Cortex" = "V1C",
  "Occipital.Pole" = "OL",
  "Lateral.Occipital.Cortex.superior.division" = "OL",
  "Lateral.Occipital.Cortex.inferior.division" = "OL",
  "Occipital.Fusiform.Gyrus" = "OL",
  "Temporal.Occipital.Fusiform.Cortex" = "OL",
  
  # Temporal
  "Inferior.Temporal.Gyrus.anterior.division" = "ITC",
  "Inferior.Temporal.Gyrus.posterior.division" = "ITC",
  "Inferior.Temporal.Gyrus.temporooccipital.part" = "ITC",
  "Temporal.Fusiform.Cortex.anterior.division" = "ITC",
  "Temporal.Fusiform.Cortex.posterior.division" = "ITC",
  "Superior.Temporal.Gyrus.anterior.division" = "STC",
  "Superior.Temporal.Gyrus.posterior.division" = "STC",
  "Middle.Temporal.Gyrus.anterior.division" = "TL",
  "Middle.Temporal.Gyrus.posterior.division" = "TL",
  "Middle.Temporal.Gyrus.temporooccipital.part" = "TL",
  "Planum.Polare" = "A1C",
  "Planum.Temporale" = "A1C",
  "Heschl.s.Gyrus.includes.H1.and.H2" = "A1C",
  
  # Frontal
  "Frontal.Orbital.Cortex" = "OFC",
  "Frontal.Operculum.Cortex" = "VFC",
  "Inferior.Frontal.Gyrus.pars.opercularis" = "VFC",
  "Inferior.Frontal.Gyrus.pars.triangularis" = "VFC",
  "Central.Opercular.Cortex" = "PL",
  "Frontal.Medial.Cortex" = "MFC",
  "Paracingulate.Gyrus" = "MFC",
  "Subcallosal.Cortex" = "MFC",
  "Cingulate.Gyrus.anterior.division" = "CNG",
  "Cingulate.Gyrus.posterior.division" = "CNG",
  "Frontal.Pole" = "FL",
  "Superior.Frontal.Gyrus" = "DFC",
  "Middle.Frontal.Gyrus" = "DFC",
  
  # Insula
  "Insular.Cortex" = "S1C"
)

# Annotate regions with macro-labels
data <- data %>%
  mutate(Region_cerebroViz = region_map[Region])

# Save annotated table
write.table(data, "FAST_plot_region_anno_v2.txt", quote = FALSE, sep = '\t', row.names = FALSE)

# Create output folder for plots
dir.create("FAST_plot_V2", showWarnings = FALSE)

# Helper for repeated cerebroViz calls
generate_cerebro_plot <- function(side_val, predictor_val, suffix) {
  dt_plot <- data[data$side == side_val & data$Predictor == predictor_val, ]
  dt_plot <- dt_plot[!is.na(dt_plot$p_value) & !is.na(dt_plot$Region_cerebroViz), ]
  dt_plot <- dt_plot[order(dt_plot$p_value), ]
  dt_plot <- dt_plot[!duplicated(dt_plot$Region_cerebroViz), c("Region_cerebroViz", "t_value")]
  
  mat_plot <- as.matrix(dt_plot[, "t_value"])
  rownames(mat_plot) <- dt_plot$Region_cerebroViz
  colnames(mat_plot) <- predictor_val
  attributes(mat_plot)$class <- "matrix"
  
  # Add dummy values for fixed color scale
  mat_plot_with_range <- rbind(
    mat_plot,
    .min_dummy = -4,
    .max_dummy = 4
  )
  attributes(mat_plot_with_range)$class <- "matrix"
  
  # Plot
  cerebroViz(mat_plot_with_range,
             filePrefix = paste0("FAST_plot_V2/FAST_", side_val, "_", suffix, "_"),
             timePoint = 1,
             naHatch = FALSE,
             regLabel = TRUE,
             palette = c("#6ff7aa", "white", "#a357e6"),
             secPalette = c("grey80", "black", "white"))
}

# ---- Generate 4 plots separately with full control ----

generate_cerebro_plot("right", "score_leading", "leading")
generate_cerebro_plot("left", "score_leading", "leading")
generate_cerebro_plot("left", "score_alternative", "alternative")
generate_cerebro_plot("right", "score_alternative", "alternative")
