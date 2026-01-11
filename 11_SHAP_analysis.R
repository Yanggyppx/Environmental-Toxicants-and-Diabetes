# 11_SHAP_analysis.R
library(shapviz)
library(fastshap)
library(ggplot2)
library(xgboost)

#-------------------------------------------------------------------------------
# 1. Load Model and Data
#-------------------------------------------------------------------------------

model <- readRDS("results/ML_models.rds")
Test_set <- readRDS("data/Test_set.rds")

# Extract XGBoost model
xgb_model <- model$XGBoost

#-------------------------------------------------------------------------------
# 2. Prepare Data
#-------------------------------------------------------------------------------

shap_data <- as.data.frame(Test_set)

if (!is.null(xgb_model$feature_names)) {
  shap_data <- shap_data[, xgb_model$feature_names, drop = FALSE]
}

cat("Samples:", nrow(shap_data), "\n")
cat("Features:", ncol(shap_data), "\n")

#-------------------------------------------------------------------------------
# 3. Define Prediction Wrapper
#-------------------------------------------------------------------------------

pred_wrapper <- function(model, newdata) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (!is.null(model$feature_names)) {
    newdata <- newdata[, model$feature_names, drop = FALSE]
  }
  as.numeric(predict(model, newdata = data.matrix(newdata)))
}

#-------------------------------------------------------------------------------
# 4. Calculate SHAP Values
#-------------------------------------------------------------------------------

cat("\nCalculating SHAP values (this may take a few minutes)...\n")

set.seed(123)

shap_values <- explain(
  object = xgb_model,
  X = shap_data,
  pred_wrapper = pred_wrapper,
  nsim = 100,
  adjust = TRUE
)

# Create shapviz object
baseline <- mean(pred_wrapper(xgb_model, shap_data))
shp <- shapviz(object = shap_values, X = shap_data, baseline = baseline)

#-------------------------------------------------------------------------------
# 5. Visualize SHAP Results
#-------------------------------------------------------------------------------

# Beeswarm plot
p_beeswarm <- sv_importance(shp, kind = "beeswarm", max_display = 15) +
  labs(
    title = "SHAP Summary Plot",
    x = "SHAP value (impact on model output)",
    y = "Features"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold")
  ) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Feature\nvalue")

ggsave("figures/SHAP_beeswarm.pdf", p_beeswarm, width = 8, height = 6)

# Bar plot
p_bar <- sv_importance(shp, kind = "bar", max_display = 15) +
  labs(
    title = "SHAP Feature Importance",
    x = "Mean |SHAP value|",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

ggsave("figures/SHAP_importance.pdf", p_bar, width = 6, height = 5)

cat("\nSHAP analysis completed!\n")