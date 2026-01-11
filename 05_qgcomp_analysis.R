#05_qgcomp_analysis.R
library(qgcomp)
library(dplyr)
library(ggplot2)

#-------------------------------------------------------------------------------
# 1. Load and Prepare Data
#-------------------------------------------------------------------------------

data <- readRDS("data/cleaned_data.rds")

# Ensure outcome is factor
data$Diabetes <- factor(data$Diabetes)

#-------------------------------------------------------------------------------
# 2. Standardize Exposure Variables
#-------------------------------------------------------------------------------

# Prioritized exposures
exposure_vars <- c("EthyleneOxide", "NAC_cHPM", "Glycidamide", "U", "Sb")

# Create standardized versions
data_scaled <- data %>%
  mutate(
    z_EthyleneOxide = scale(EthyleneOxide),
    z_NAC_cHPM = scale(NAC_cHPM),
    z_Glycidamide = scale(Glycidamide),
    z_U = scale(U),
    z_Sb = scale(Sb)
  )

# Exposure names for qgcomp
exp_names <- c("z_EthyleneOxide", "z_NAC_cHPM", "z_Glycidamide", "z_U", "z_Sb")

#-------------------------------------------------------------------------------
# 3. Run qgcomp Without Bootstrap (Quick Estimation)
#-------------------------------------------------------------------------------

cat("===== Running qgcomp (non-bootstrap) =====\n")

result_noboot <- qgcomp.noboot(
  Diabetes ~ z_EthyleneOxide + z_NAC_cHPM + z_Glycidamide + z_U + z_Sb +
    Age + BMI + Sex + Education + MaritalStatus + PIR + Alcohol,
  expnms = exp_names,
  data = data_scaled,
  family = binomial(),
  q = 4
)

cat("\nNon-bootstrap Results:\n")
print(summary(result_noboot))

# Extract weights
weights_noboot <- data.frame(
  Exposure = names(result_noboot$pos.weights),
  Weight = as.numeric(result_noboot$pos.weights),
  Direction = "Positive"
)

# Add negative weights if any
if (length(result_noboot$neg.weights) > 0) {
  neg_weights <- data.frame(
    Exposure = names(result_noboot$neg.weights),
    Weight = as.numeric(result_noboot$neg.weights),
    Direction = "Negative"
  )
  weights_noboot <- rbind(weights_noboot, neg_weights)
}

#-------------------------------------------------------------------------------
# 4. Run qgcomp With Bootstrap (For Inference)
#-------------------------------------------------------------------------------

cat("\n===== Running qgcomp (bootstrap, B=1000) =====\n")
cat("This may take several minutes...\n")

set.seed(125)

result_boot <- qgcomp.boot(
  Diabetes ~ z_EthyleneOxide + z_NAC_cHPM + z_Glycidamide + z_U + z_Sb +
    Age + BMI + Sex + Education + MaritalStatus + PIR + Alcohol,
  expnms = exp_names,
  data = data_scaled,
  family = binomial(),
  q = 4,
  B = 1000,
  seed = 125,
  rr = FALSE  # Return OR instead of RR
)

cat("\nBootstrap Results:\n")
print(summary(result_boot))

#-------------------------------------------------------------------------------
# 5. Extract and Format Results
#-------------------------------------------------------------------------------

# Mixture effect (psi1)
mixture_effect <- data.frame(
  OR = exp(result_boot$coef["psi1"]),
  Lower_CI = exp(result_boot$ci[1]),
  Upper_CI = exp(result_boot$ci[2]),
  P_value = result_boot$pval["psi1"]
)

cat("\nMixture Effect (per quartile increase):\n")
cat(
  "OR:", round(mixture_effect$OR, 3),
  "(95% CI:", round(mixture_effect$Lower_CI, 3), "-",
  round(mixture_effect$Upper_CI, 3), ")\n"
)

#-------------------------------------------------------------------------------
# 6. Visualize Weights
#-------------------------------------------------------------------------------

# Prepare data for plotting
weights_noboot$Exposure <- gsub("^z_", "", weights_noboot$Exposure)
weights_noboot <- weights_noboot %>% arrange(Weight)
weights_noboot$Exposure <- factor(
  weights_noboot$Exposure,
  levels = weights_noboot$Exposure
)

p_weights <- ggplot(weights_noboot, aes(x = Weight, y = Exposure)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_col(aes(fill = Direction), width = 0.6) +
  scale_fill_manual(
    values = c("Positive" = "#B2182B", "Negative" = "#2166AC")
  ) +
  labs(
    title = "qgcomp Exposure Weights",
    x = "Weight",
    y = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

ggsave("figures/qgcomp_weights.pdf",
       p_weights, width = 6, height = 4)

#-------------------------------------------------------------------------------
# 7. Default qgcomp Plot
#-------------------------------------------------------------------------------

pdf("figures/qgcomp_mixture_plot.pdf", width = 8, height = 6)
plot(result_boot)
dev.off()

#-------------------------------------------------------------------------------
# 8. Subgroup Analysis (Optional)
#-------------------------------------------------------------------------------

cat("\n===== Subgroup Analysis =====\n")

# Define subgroups
subgroups <- list(
  BMI_high = data_scaled %>% filter(BMI >= 25),
  BMI_low  = data_scaled %>% filter(BMI < 25),
  Age_old  = data_scaled %>% filter(Age >= 60),
  Age_young = data_scaled %>% filter(Age < 60)
)

subgroup_results <- list()

for (subgroup_name in names(subgroups)) {

  cat("Analyzing subgroup:", subgroup_name, "... ")

  subgroup_data <- subgroups[[subgroup_name]]

  # Skip if too few samples
  if (nrow(subgroup_data) < 100) {
    cat("Skipped (n < 100)\n")
    next
  }

  # Remove BMI or Age from covariates if stratified
  if (grepl("BMI", subgroup_name)) {
    formula <- Diabetes ~ z_EthyleneOxide + z_NAC_cHPM +
      z_Glycidamide + z_U + z_Sb +
      Age + Sex + Education + MaritalStatus + PIR + Alcohol
  } else if (grepl("Age", subgroup_name)) {
    formula <- Diabetes ~ z_EthyleneOxide + z_NAC_cHPM +
      z_Glycidamide + z_U + z_Sb +
      BMI + Sex + Education + MaritalStatus + PIR + Alcohol
  }

  tryCatch({
    res <- qgcomp.boot(
      formula,
      expnms = exp_names,
      data = subgroup_data,
      family = binomial(),
      q = 4,
      B = 500,
      seed = 125,
      rr = FALSE
    )

    subgroup_results[[subgroup_name]] <- data.frame(
      Subgroup = subgroup_name,
      N = nrow(subgroup_data),
      OR = exp(res$coef["psi1"]),
      Lower_CI = exp(res$ci[1]),
      Upper_CI = exp(res$ci[2])
    )

    cat("OR =", round(exp(res$coef["psi1"]), 3), "\n")

  }, error = function(e) {
    cat("Failed\n")
  })
}

# Combine subgroup results
if (length(subgroup_results) > 0) {
  subgroup_table <- do.call(rbind, subgroup_results)
  write.csv(
    subgroup_table,
    "results/qgcomp_subgroup_results.csv",
    row.names = FALSE
  )
  cat("\nSubgroup results saved.\n")
}

#-------------------------------------------------------------------------------
# 9. Save All Results
#-------------------------------------------------------------------------------

results_summary <- list(
  noboot = result_noboot,
  boot = result_boot,
  weights = weights_noboot,
  mixture_effect = mixture_effect
)

saveRDS(results_summary, "results/qgcomp_results.rds")
write.csv(weights_noboot, "results/qgcomp_weights.csv", row.names = FALSE)
write.csv(mixture_effect, "results/qgcomp_mixture_effect.csv", row.names = FALSE)

cat("\n===== qgcomp Analysis Completed =====\n")