# ==============================================================================
# 1.Fig3A.R
# Author: Xiaopu Zhou
# Purpose:
#   Analyze APOE expression effects across ancestry groups:
#     - Visualize ancestry-stratified Delta vs APOE effect
#     - Fit Bayesian measurement error (Deming) regression
#     - Report Bayesian p-values and correlation
# ==============================================================================



# ==== [1] Load Libraries ====
library(ggplot2)
library(dplyr)
library(MethComp)
library(brms)
library(posterior)  # for as_draws_df()

# ==== [2] Define Color Scheme for Ancestries ====
ancestry_colors <- c(
  AFR = "#B22222",  # Firebrick
  AMR = "#1E88A8",  # Blue
  EAS = "#00775E",  # Teal
  EUR = "#2B3E68",  # Indigo
  SAS = "#C46245"   # Burnt sienna
)

# ==== [3] Define Input Data ====

# Delta scores (evolution score)
x_data <- read.table(text = "
Superpopulation_code mean_delta sem_delta
AFR -0.00038 0.0000301
AMR -0.000416 0.0000588
EAS -0.000223 0.0000181
EUR -0.00027 0.0000328
SAS -0.000308 0.0000289", header = TRUE)

# APOE expression effect size
y_data <- read.table(text = "
Superpopulation_code mean_APOE sem_APOE
AFR -0.3210744 0.06671819
AMR -0.2959885 0.1241098
EAS  0.2260632 0.08085577
EUR  0.2129778 0.09709001
SAS  0.1532142 0.10502246", header = TRUE)

# Merge input for plotting and modeling
plot_data <- left_join(x_data, y_data, by = "Superpopulation_code")

# ==== [4] Plot Ancestry-Stratified Correlation with Error Bars ====
p <- ggplot(plot_data, aes(x = mean_delta, y = mean_APOE, color = Superpopulation_code)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_APOE - sem_APOE, ymax = mean_APOE + sem_APOE),
                width = 0.00001, size = 0.8) +
  geom_errorbarh(aes(xmin = mean_delta - sem_delta, xmax = mean_delta + sem_delta),
                 height = 0.05, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "#D3D3D3",
              linetype = "dashed", size = 0.6, alpha = 0.3) +
  scale_color_manual(values = ancestry_colors) +
  theme_classic(base_size = 12) +
  labs(x = "Mean Delta (± SE)", y = "Mean APOE Effect (± SE)") +
  theme(legend.title = element_blank())

# Save plot to PNG
png("APOE_ancestry_expression_correlation.png", res = 1200,
    width = 12, height = 9, pointsize = 10, units = "cm")
print(p)
dev.off()

# ==== [5] Prepare Data for Bayesian Measurement Error Regression ====
merged_data <- merge(x_data, y_data, by = "Superpopulation_code")
bayes_data <- data.frame(
  x = merged_data$mean_delta,
  y = merged_data$mean_APOE,
  x_se = merged_data$sem_delta,
  y_se = merged_data$sem_APOE
)

# ==== [6] Fit Bayesian Measurement Error Model (Deming Regression) ====
fit_bayes <- brm(
  bf(y | se(y_se) ~ me(x, x_se)),
  data = bayes_data,
  family = gaussian(),
  chains = 4,
  iter = 6000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 123
)

print(summary(fit_bayes))

# ==== [7] Compute Posterior Probability (Bayesian p-value) ====
posterior_samples <- as_draws_df(fit_bayes)
slope_col <- grep("^bsp_", names(posterior_samples), value = TRUE)

p_gt_0 <- mean(posterior_samples[[slope_col]] > 0)
p_lt_0 <- mean(posterior_samples[[slope_col]] < 0)
bayes_p <- 2 * min(p_gt_0, p_lt_0)

cat("P(beta > 0):", round(p_gt_0, 4), "\n")
cat("P(beta < 0):", round(p_lt_0, 4), "\n")
cat("Two-tailed Bayesian p-value:", round(bayes_p, 6), "\n")

# ==== [8] Compute Bayesian Correlation Coefficient (R) ====
fitted_matrix <- fitted(fit_bayes, summary = FALSE)
y_obs <- fit_bayes$data$y

R_vals <- apply(fitted_matrix, 1, function(yhat) cor(yhat, y_obs))
R_estimate <- mean(R_vals)
R_ci <- quantile(R_vals, probs = c(0.025, 0.975))

cat("Bayesian R estimate:", round(R_estimate, 4), "\n")
cat("95% Credible Interval:", round(R_ci[1], 4), "-", round(R_ci[2], 4), "\n")
