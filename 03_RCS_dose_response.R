# 03_RCS_dose_response.R
library(rms)
library(ggplot2)
library(gridExtra)

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

data <- readRDS("data/cleaned_data.rds")

# Set up data distribution for rms package
dd <- datadist(data)
options(datadist = "dd")

#-------------------------------------------------------------------------------
# 2. Define Exposures and Covariates
#-------------------------------------------------------------------------------

# Prioritized exposures from DSA analysis
exposures <- c("EthyleneOxide", "NAC_cHPM", "Glycidamide", "U", "Sb")

# Covariates for adjustment
covariates <- "Age + BMI + Sex + Education + MaritalStatus + PIR + Alcohol"

#-------------------------------------------------------------------------------
# 3. Function for RCS Analysis and Visualization
#-------------------------------------------------------------------------------

run_rcs_analysis <- function(data, exposure, outcome = "Diabetes",
                             covariates, knots = 3) {

  # Build formula
  formula_str <- paste0(
    outcome, " ~ rcs(", exposure, ", ", knots, ") + ", covariates
  )
  formula <- as.formula(formula_str)

  # Fit logistic regression model
  fit <- lrm(formula, data = data)

  # ANOVA for significance testing
  an <- anova(fit)

  # Extract p-values
  total_p <- an[which(rownames(an) == exposure), "P"]
  nonlinear_p <- an[which(rownames(an) == " Nonlinear"), "P"]

  # Generate predictions
  OR_pred <- Predict(fit, var = exposure, fun = exp, ref.zero = TRUE)
  OR_pred <- as.data.frame(OR_pred)
  colnames(OR_pred)[1] <- "exposure_value"

  # Create plot
  p <- ggplot() +
    geom_line(
      data = OR_pred,
      aes(x = exposure_value, y = yhat),
      linetype = 1, linewidth = 0.5,
      alpha = 1, colour = "#4692C5"
    ) +
    geom_ribbon(
      data = OR_pred,
      aes(x = exposure_value, ymin = lower, ymax = upper),
      alpha = 0.2, fill = "#4692C5"
    ) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste("P for overall =", format(total_p, digits = 3)),
      hjust = 1.1, vjust = 1.5, family = "serif", size = 3
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste("P for nonlinear =", format(nonlinear_p, digits = 3)),
      hjust = 1.1, vjust = 3, family = "serif", size = 3
    ) +
    theme_classic() +
    labs(x = "Concentration", y = "OR (95% CI)") +
    ggtitle(exposure) +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(
        hjust = 0.5, size = 10, face = "bold"
      ),
      axis.line = element_line(linewidth = 0.3)
    )

  # Return results
  return(list(
    model = fit,
    anova = an,
    total_p = total_p,
    nonlinear_p = nonlinear_p,
    plot = p
  ))
}

#-------------------------------------------------------------------------------
# 4. Run RCS Analysis for All Exposures
#-------------------------------------------------------------------------------

cat("===== Running RCS Analysis =====\n\n")

rcs_results <- list()
plot_list <- list()

for (exp in exposures) {

  cat("Analyzing:", exp, "... ")

  tryCatch({

    result <- run_rcs_analysis(
      data = data,
      exposure = exp,
      covariates = covariates
    )

    rcs_results[[exp]] <- result
    plot_list[[exp]] <- result$plot

    cat(
      "P_overall =", format(result$total_p, digits = 3),
      ", P_nonlinear =", format(result$nonlinear_p, digits = 3), "\n"
    )

  }, error = function(e) {
    cat("Failed:", conditionMessage(e), "\n")
  })
}

#-------------------------------------------------------------------------------
# 5. Combine and Save Plots
#-------------------------------------------------------------------------------

combined_plot <- do.call(
  grid.arrange,
  c(plot_list, ncol = 3)
)

ggsave(
  "figures/RCS_dose_response_all.pdf",
  combined_plot,
  width = 12, height = 8, dpi = 300
)

#-------------------------------------------------------------------------------
# 6. Export Summary Table
#-------------------------------------------------------------------------------

summary_table <- data.frame(
  Exposure = names(rcs_results),
  P_overall = sapply(rcs_results, function(x) x$total_p),
  P_nonlinear = sapply(rcs_results, function(x) x$nonlinear_p),
  Significant = sapply(
    rcs_results,
    function(x) ifelse(x$total_p < 0.05, "Yes", "No")
  ),
  Nonlinear = sapply(
    rcs_results,
    function(x) ifelse(x$nonlinear_p < 0.05, "Yes", "No")
  )
)

write.csv(
  summary_table,
  "results/RCS_summary.csv",
  row.names = FALSE
)

cat("\n===== RCS Analysis Completed =====\n")
print(summary_table)
