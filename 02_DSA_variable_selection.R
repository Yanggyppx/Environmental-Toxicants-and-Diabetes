#02_DSA_variable_selection.R
library(DSA)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

data <- readRDS("data/cleaned_data.rds")

#-------------------------------------------------------------------------------
# 2. Define Variables
#-------------------------------------------------------------------------------

# Outcome variable
Y <- "Diabetes"

# Fixed covariates (confounders)
fixed_covariates <- c(
  "Age", "BMI", "Sex", "Race",
  "Education", "MaritalStatus", "PIR",
  "Smoking", "Alcohol"
)

# Exposure variables (all other variables)
exposure_vars <- setdiff(colnames(data), c(Y, fixed_covariates))

cat("Number of exposure variables:", length(exposure_vars), "\n")

#-------------------------------------------------------------------------------
# 3. Run DSA with 100 Iterations
#-------------------------------------------------------------------------------

base_formula <- as.formula(paste(Y, "~ 1"))
selected_vars <- list()

set.seed(123)
n_repeats <- 100

for (i in 1:n_repeats) {

  cat("Iteration:", i, "/", n_repeats, "\r")

  res <- DSA(
    formula = base_formula,
    data = data[, c(Y, fixed_covariates, exposure_vars)],
    family = binomial(),
    force.in = fixed_covariates,
    maxsize = 15,
    maxorderint = 1,
    maxsumofpow = 1,
    cv = TRUE,
    nfolds = 5
  )

  # Extract selected variables
  selected <- attr(terms(res$model.selected), "term.labels")
  clean_vars <- gsub("^I\\((.+)\\^\\d+\\)$", "\\1", selected)
  clean_vars <- gsub("\\d+$", "", clean_vars)
  selected_exposures <- setdiff(clean_vars, fixed_covariates)
  selected_vars[[i]] <- selected_exposures
}

cat("\n")

#-------------------------------------------------------------------------------
# 4. Summarize Selection Frequency
#-------------------------------------------------------------------------------

selected_flat <- unlist(selected_vars)
selection_freq <- sort(
  table(selected_flat) / n_repeats,
  decreasing = TRUE
)

# Filter variables selected in >= 5% of iterations
final_exposures <- names(selection_freq[selection_freq >= 0.05])

cat("\nVariables selected in >= 5% of iterations:\n")
print(selection_freq[selection_freq >= 0.05])

#-------------------------------------------------------------------------------
# 5. Save Results
#-------------------------------------------------------------------------------

results <- list(
  all_selections = selected_vars,
  selection_frequency = selection_freq,
  final_exposures = final_exposures
)

saveRDS(results, "results/DSA_results.rds")
write.csv(
  as.data.frame(selection_freq),
  "results/DSA_selection_frequency.csv"
)

cat("\nDSA analysis completed!\n")
cat(
  "Prioritized exposures:",
  paste(final_exposures, collapse = ", "),
  "\n"
)
