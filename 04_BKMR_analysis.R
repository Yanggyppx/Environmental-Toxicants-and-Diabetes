#04_BKMR_analysis.R
library(bkmr)
library(ggplot2)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Load and Prepare Data
#-------------------------------------------------------------------------------

data <- readRDS("data/cleaned_data.rds")

# Define exposure and covariate matrices
# Adjust column indices based on your prioritized exposures
exposure_cols <- c("EthyleneOxide", "NAC_cHPM", "Glycidamide", "U", "Sb")
covar_cols <- c("Age", "BMI", "Sex", "Race", "Education",
                "MaritalStatus", "PIR", "Smoking", "Alcohol")

expos <- data.matrix(data[, exposure_cols])
covar <- data.matrix(data[, covar_cols])
Y <- data$Diabetes

# Standardize exposures
scale_expos <- scale(expos)

#-------------------------------------------------------------------------------
# 2. Fit BKMR Model
#-------------------------------------------------------------------------------

set.seed(1234)

# Use knots for computational efficiency with large samples
knots50 <- fields::cover.design(scale_expos, nd = 50)$design

fit_bkmr <- kmbayes(
  y = Y,
  Z = scale_expos,
  X = covar,
  iter = 10000,
  verbose = FALSE,
  varsel = TRUE,
  family = "binomial",
  est.h = TRUE,
  knots = knots50
)

#-------------------------------------------------------------------------------
# 3. Extract Posterior Inclusion Probabilities (PIPs)
#-------------------------------------------------------------------------------

PIPs <- ExtractPIPs(fit_bkmr)
cat("\nPosterior Inclusion Probabilities:\n")
print(PIPs)

write.csv(PIPs, "results/BKMR_PIPs.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# 4. Univariate Exposure-Response
#-------------------------------------------------------------------------------

pred_univar <- PredictorResponseUnivar(fit = fit_bkmr)

p1 <- ggplot(pred_univar,
             aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ variable) +
  ylab("h(z)") +
  theme_bw() +
  theme(text = element_text(family = "serif"))

ggsave("figures/BKMR_univariate_response.pdf",
       p1, width = 8, height = 6)

#-------------------------------------------------------------------------------
# 5. Overall Mixture Effect
#-------------------------------------------------------------------------------

risks_overall <- OverallRiskSummaries(
  fit = fit_bkmr,
  qs = seq(0.25, 0.75, by = 0.05),
  q.fixed = 0.5
)

p2 <- ggplot(risks_overall,
             aes(quantile, est,
                 ymin = est - 1.96 * sd,
                 ymax = est + 1.96 * sd)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  xlab("Quantile") +
  ylab("Estimate") +
  theme_bw() +
  theme(
    text = element_text(family = "serif", size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )

ggsave("figures/BKMR_overall_effect.pdf",
       p2, width = 6, height = 4)

#-------------------------------------------------------------------------------
# 6. Save Model
#-------------------------------------------------------------------------------

saveRDS(fit_bkmr, "results/BKMR_model.rds")

cat("\nBKMR analysis completed!\n")
