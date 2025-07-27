# ==============================================================================
# 6.Fig3F.R
# Author: Xiaopu Zhou
# Purpose:
#   - Merges normalized EVO2 score and APOE-ε4 odds ratio across ancestries
#   - Visualizes their correlation with uncertainty bars
#   - Performs Bayesian Deming regression with brms
#   - Outputs Bayesian p-value and correlation
# ==============================================================================

# ==== [1] Set working directory ====
setwd("C:/Users/Fred/OneDrive/2025April/EVO2_APOE/Figure/Re_Figure/Figure3_Pangenome/panel_f")

# ==== [2] Load libraries ====
library(ggplot2)
library(dplyr)
library(brms)
library(posterior)

# ==== [3] Define ancestry-specific color palette ====
npg_color_lines <- c(
  AFR = "#B22222",  # Firebrick
  AMR = "#1E88A8",  # Blue
  EAS = "#00775E",  # Teal
  EUR = "#2B3E68",  # Indigo
  SAS = "#C46245"   # Burnt sienna
)

# ==== [4] Load input data: mean and SE of EVO2 scores and ORs ====

# Normalized EVO2 score across populations
x_data <- read.table(text = "
Superpopulation_code mean_EVO2 sem_EVO2
SAS -0.000124 0.00007
EUR 0.000049 0.0000653
EAS 0.000133 0.0000299
AMR -0.000229 0.000112
AFR -0.000252 0.0000605", header = TRUE)

# APOE-ε4 odds ratios (from external association results)
y_data <- read.table(text = "
Superpopulation_code mean_OR sem_OR
EUR 3.46 0.096938776
EAS 4.54 0.301020408
AMR 1.9 0.135204082
AFR 2.18 0.150510204", header = TRUE)

# ==== [5] Merge data and filter for valid OR ====
plot_data <- left_join(x_data, y_data, by = "Superpopulation_code")

# Drop SAS (no OR value available)
bayes_data <- plot_data[complete.cases(plot_data), ]

# ==== [6] Visualization: Scatter plot with error bars and regression ====
p <- ggplot(plot_data, aes(x = mean_EVO2, y = mean_OR, color = Superpopulation_code)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_OR - sem_OR, ymax = mean_OR + sem_OR),
                width = 0.00001, size = 0.8) +
  geom_errorbarh(aes(xmin = mean_EVO2 - sem_EVO2, xmax = mean_EVO2 + sem_EVO2),
                 height = 0.1, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "#D3D3D3",
              linetype = "dashed", size = 0.6, alpha = 0.3) +
  scale_color_manual(values = npg_color_lines) +
  theme_classic(base_size = 12) +
  labs(x = "Normalized EVO2 Score (± SE)", y = "APOE-ε4 Odds Ratio (± SE)") +
  theme(legend.title = element_blank())

# Save the figure
png(file = "EVO2_OR_correlation.png", res = 1200, width = 12, height = 9, pointsize = 10, units = "cm")
print(p)
dev.off()

# ==== [7] Prepare input for Bayesian Deming regression ====
bayes_input <- data.frame(
  x = bayes_data$mean_EVO2,
  y = bayes_data$mean_OR,
  x_se = bayes_data$sem_EVO2,
  y_se = bayes_data$sem_OR
)

# ==== [8] Fit Bayesian model using brms ====
fit_bayes <- brm(
  bf(y | se(y_se) ~ me(x, x_se)),
  data = bayes_input,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 123
)

# Print model summary
summary(fit_bayes)

# ==== [9] Extract posterior draws and test slope directionality ====
posterior_draws <- as_draws_df(fit_bayes)
slope_col <- grep("^bsp_", names(posterior_draws), value = TRUE)

p_gt_0 <- mean(posterior_draws[[slope_col]] > 0)
p_lt_0 <- mean(posterior_draws[[slope_col]] < 0)
bayes_p <- 2 * min(p_gt_0, p_lt_0)

cat("P(beta > 0):", round(p_gt_0, 4), "\n")
cat("P(beta < 0):", round(p_lt_0, 4), "\n")
cat("Two-tailed Bayesian p-value:", round(bayes_p, 6), "\n")

# ==== [10] Estimate Bayesian correlation coefficient ====
fitted_matrix <- fitted(fit_bayes, summary = FALSE)
y_obs <- fit_bayes$data$y

R_vals <- apply(fitted_matrix, 1, function(yhat) cor(yhat, y_obs))
R_estimate <- mean(R_vals)
R_ci <- quantile(R_vals, c(0.025, 0.975))

cat("Bayesian R estimate:", round(R_estimate, 4), "\n")
cat("95% Credible Interval:", round(R_ci[1], 4), "-", round(R_ci[2], 4), "\n")
                  