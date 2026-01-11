# 09_ML_model_training.R
# Load required packages
library(caret)
library(glmnet)
library(randomForestSRC)
library(xgboost)
library(pROC)
library(dplyr)

# Load custom ML functions
source("05_Machine_Learning/refer.ML.R")

#-------------------------------------------------------------------------------
# 1. Load and Prepare Data
#-------------------------------------------------------------------------------

cat("===== Loading Data =====\n")

train_data <- read.csv("data/train_data.csv")
test_data <- read.csv("data/test_data.csv")

# Separate features and labels
Train_expr <- train_data[, 3:ncol(train_data)]
Train_class <- factor(train_data$Group)

Test_expr <- test_data[, 3:ncol(test_data)]
Test_class <- factor(test_data$Group)

cat("Training set:", nrow(Train_expr), "samples x", ncol(Train_expr), "features\n")
cat("Test set:", nrow(Test_expr), "samples x", ncol(Test_expr), "features\n")

#-------------------------------------------------------------------------------
# 2. Data Preprocessing
#-------------------------------------------------------------------------------

cat("\n===== Preprocessing =====\n")

# Align features
common_features <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[, common_features])
Test_expr <- as.matrix(Test_expr[, common_features])

# Remove zero-variance features
feature_sd <- apply(Train_expr, 2, sd)
valid_features <- feature_sd > 0
Train_expr <- Train_expr[, valid_features]
Test_expr <- Test_expr[, valid_features]

# Standardize
train_mean <- apply(Train_expr, 2, mean)
train_sd <- apply(Train_expr, 2, sd)

Train_set <- scale(Train_expr, center = train_mean, scale = train_sd)
Test_set <- scale(Test_expr, center = train_mean, scale = train_sd)

cat("Final features:", ncol(Train_set), "\n")

#-------------------------------------------------------------------------------
# 3. Feature Pre-selection (Prevent Overfitting)
#-------------------------------------------------------------------------------

cat("\n===== Feature Selection =====\n")

n_samples <- nrow(Train_set)
max_features <- max(15, floor(n_samples / 5))

if (ncol(Train_set) > max_features) {
  
  # Univariate significance screening
  Train_class_num <- as.numeric(Train_class) - 1
  
  pvals <- apply(Train_set, 2, function(x) {
    tryCatch(wilcox.test(x ~ Train_class_num)$p.value, error = function(e) 1)
  })
  
  selected_features <- names(sort(pvals)[1:max_features])
  
  Train_set <- Train_set[, selected_features]
  Test_set <- Test_set[, selected_features]
  
  cat("Selected", length(selected_features), "features\n")
}

#-------------------------------------------------------------------------------
# 4. Train Models
#-------------------------------------------------------------------------------

cat("\n===== Training Models =====\n")

# Load method list
methods <- read.table("05_Machine_Learning/refer.methodLists.txt", 
                      header = TRUE, sep = "\t")$Model
methods <- gsub("-| ", "", methods)

# Prepare labels
Train_label <- data.frame(Type = Train_class)
Train_label_num <- data.frame(Type = as.numeric(Train_class) - 1)

# Train models
model <- list()
set.seed(123)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  cat(sprintf("[%d/%d] %s... ", i, length(methods), method))
  
  tryCatch({
    
    method_split <- strsplit(method, "\\+")[[1]]
    if (length(method_split) == 1) method_split <- c("simple", method_split)
    
    # Select appropriate label format
    use_numeric <- grepl("GBM|XGBoost|plsRglm", method_split[2])
    current_label <- if (use_numeric) Train_label_num else Train_label
    
    # Train model
    model[[method]] <- RunML(
      method = method_split[2],
      Train_set = Train_set,
      Train_label = current_label,
      mode = "Model",
      classVar = "Type"
    )
    
    cat("Success\n")
    
  }, error = function(e) {
    cat("Failed\n")
  })
}

cat("\nSuccessfully trained:", length(model), "models\n")

#-------------------------------------------------------------------------------
# 5. Save Models
#-------------------------------------------------------------------------------

saveRDS(model, "results/ML_models.rds")

cat("\nModel training completed!\n")