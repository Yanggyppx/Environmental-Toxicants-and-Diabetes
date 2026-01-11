#01_ExWAS_analysis.R
library(dplyr)
library(broom)
library(purrr)
library(readxl)

# Step 1: Read data
data <- read_excel("data/exposure_diabetes_log_transformed.xlsx")

# Step 2: Convert variable types
data[, 18:63] <- lapply(data[, 18:63], function(x) as.numeric(as.character(x)))
data[[8]] <- as.numeric(as.character(data[[8]]))   # Age
data[[9]] <- as.numeric(as.character(data[[9]]))   # BMI

cat_vars <- colnames(data)[10:17]
data[, cat_vars] <- lapply(data[, cat_vars], as.factor)

# Step 3: Define variables
exposure_vars <- colnames(data)[18:63]
outcome_var <- colnames(data)[7]
covariates <- colnames(data)[8:17]

# Step 4: ExWAS analysis
covariate_string <- paste(covariates, collapse = " + ")
error_vars <- c()
empty_result_vars <- c()

exwas_result <- map_dfr(exposure_vars, function(exp) {
  
  formula <- as.formula(
    paste0(outcome_var, " ~ ", exp, " + ", covariate_string)
  )
  
  fit <- tryCatch(
    glm(formula, data = data, family = binomial()),
    error = function(e) {
      error_vars <<- c(error_vars, exp)
      return(NULL)
    }
  )
  
  # Model fitting failed or did not converge
  if (is.null(fit) || !fit$converged) {
    error_vars <<- c(error_vars, exp)
    return(NULL)
  }
  
  tidy_fit <- tidy(fit)
  
  # Exact term matching to ensure correct extraction of the main effect
  beta_row <- tidy_fit %>% filter(term == exp)
  if (nrow(beta_row) == 0) {
    empty_result_vars <<- c(empty_result_vars, exp)
    return(NULL)
  }
  
  beta <- beta_row$estimate
  se <- beta_row$std.error
  pval <- beta_row$p.value
  OR <- exp(beta)
  OR_lower <- exp(beta - 1.96 * se)
  OR_upper <- exp(beta + 1.96 * se)
  
  tibble(
    exposure = exp,
    beta = beta,
    se = se,
    pval = pval,
    OR = OR,
    OR_ci_lower = OR_lower,
    OR_ci_upper = OR_upper
  )
})

# Step 5: Report problematic variables
cat("Variables with model fitting errors (glm error or non-convergence):\n")
print(error_vars)

cat("Variables without extracted main effects (term matching failure):\n")
print(empty_result_vars)

# Step 6: FDR correction
exwas_result <- exwas_result %>%
  mutate(
    fdr_bh = p.adjust(pval, method = "BH"),
    sig_raw = pval < 0.05,
    sig_fdr = fdr_bh < 0.05
  )

# Step 7: Export results
write.csv(
  exwas_result,
  "results/exwas_results_unweighted.csv",
  row.names = FALSE
)