#10_model_evaluation.R
library(caret)
library(pROC)
library(ggplot2)
library(gridExtra)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Load Models and Data
#-------------------------------------------------------------------------------

cat("===== Loading Models and Data =====\n")

model <- readRDS("results/ML_models.rds")
Train_set <- readRDS("data/Train_set.rds")
Test_set <- readRDS("data/Test_set.rds")
Train_class <- readRDS("data/Train_class.rds")
Test_class <- readRDS("data/Test_class.rds")

# Load prediction results
classTab <- read.table("results/model.classMatrix.txt", 
                       header = TRUE, sep = "\t", row.names = 1)
riskTab <- read.table("results/model.riskMatrix.txt", 
                      header = TRUE, sep = "\t", row.names = 1)

n_train <- nrow(Train_set)
n_test <- nrow(Test_set)

cat("Training samples:", n_train, "\n")
cat("Test samples:", n_test, "\n")
cat("Available models:", length(model), "\n")

#-------------------------------------------------------------------------------
# 2. Select Best Model
#-------------------------------------------------------------------------------

# Read AUC matrix to find best model
auc_matrix <- read.table("results/model.AUCmatrix.txt", 
                         header = TRUE, sep = "\t", row.names = 1)
mean_auc <- rowMeans(auc_matrix)
best_model_name <- names(which.max(mean_auc))

cat("\nBest model:", best_model_name, "\n")
cat("Mean AUC:", round(max(mean_auc), 4), "\n")

#-------------------------------------------------------------------------------
# 3. Extract Predictions for Best Model
#-------------------------------------------------------------------------------

# True labels
train_true <- Train_class
test_true <- Test_class

# Predicted labels (rows 1:n_train are training, rest are test)
train_pred_class <- factor(classTab[1:n_train, best_model_name], 
                           levels = levels(train_true))
test_pred_class <- factor(classTab[(n_train+1):(n_train+n_test), best_model_name], 
                          levels = levels(test_true))

# Predicted probabilities
train_pred_prob <- riskTab[1:n_train, best_model_name]
test_pred_prob <- riskTab[(n_train+1):(n_train+n_test), best_model_name]

#-------------------------------------------------------------------------------
# 4. Confusion Matrix
#-------------------------------------------------------------------------------

cat("\n===== Confusion Matrix =====\n")

train_cm <- confusionMatrix(train_pred_class, train_true)
test_cm <- confusionMatrix(test_pred_class, test_true)

cat("\n[Training Set]\n")
print(train_cm)

cat("\n[Test Set]\n")
print(test_cm)

#-------------------------------------------------------------------------------
# 5. Confusion Matrix Visualization
#-------------------------------------------------------------------------------

plot_confusion_matrix <- function(cm, title) {
  
  cm_df <- as.data.frame(cm$table)
  
  p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = Freq), size = 6, fontface = "bold", color = "black") +
    scale_fill_gradient(low = "#E8F4F8", high = "#2E86AB") +
    labs(title = title, x = "Actual", y = "Predicted") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

p_cm_train <- plot_confusion_matrix(train_cm, "Training Set")
p_cm_test <- plot_confusion_matrix(test_cm, "Test Set")

#-------------------------------------------------------------------------------
# 6. ROC Curves
#-------------------------------------------------------------------------------

cat("\n===== ROC Analysis =====\n")

# Convert labels to numeric
train_true_num <- as.numeric(as.character(train_true))
test_true_num <- as.numeric(as.character(test_true))

# Calculate ROC
train_roc <- roc(train_true_num, train_pred_prob, quiet = TRUE)
test_roc <- roc(test_true_num, test_pred_prob, quiet = TRUE)

train_auc <- auc(train_roc)
test_auc <- auc(test_roc)

cat("Training AUC:", round(train_auc, 4), "\n")
cat("Test AUC:", round(test_auc, 4), "\n")

# ROC plot function
plot_roc <- function(roc_obj, auc_val, title, color = "#E63946") {
  
  roc_df <- data.frame(
    TPR = roc_obj$sensitivities,
    FPR = 1 - roc_obj$specificities
  )
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = color, linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray50", linewidth = 0.5) +
    annotate("text", x = 0.6, y = 0.3,
             label = paste0("AUC = ", round(auc_val, 4)),
             size = 5, fontface = "bold", color = color) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    coord_fixed() +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  
  return(p)
}

p_roc_train <- plot_roc(train_roc, train_auc, "Training Set")
p_roc_test <- plot_roc(test_roc, test_auc, "Test Set")

#-------------------------------------------------------------------------------
# 7. Combined Figure
#-------------------------------------------------------------------------------

combined_plot <- grid.arrange(
  p_cm_train, p_cm_test, p_roc_train, p_roc_test,
  ncol = 2, nrow = 2,
  top = textGrob(paste("Model Evaluation:", best_model_name),
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave("figures/Model_evaluation_combined.pdf", combined_plot, 
       width = 10, height = 10, dpi = 300)
ggsave("figures/Model_evaluation_combined.png", combined_plot, 
       width = 10, height = 10, dpi = 300)

#-------------------------------------------------------------------------------
# 8. Performance Metrics Summary
#-------------------------------------------------------------------------------

metrics_table <- data.frame(
  Dataset = c("Training", "Test"),
  N = c(n_train, n_test),
  Accuracy = c(
    round(train_cm$overall["Accuracy"], 4),
    round(test_cm$overall["Accuracy"], 4)
  ),
  Sensitivity = c(
    round(train_cm$byClass["Sensitivity"], 4),
    round(test_cm$byClass["Sensitivity"], 4)
  ),
  Specificity = c(
    round(train_cm$byClass["Specificity"], 4),
    round(test_cm$byClass["Specificity"], 4)
  ),
  PPV = c(
    round(train_cm$byClass["Pos Pred Value"], 4),
    round(test_cm$byClass["Pos Pred Value"], 4)
  ),
  NPV = c(
    round(train_cm$byClass["Neg Pred Value"], 4),
    round(test_cm$byClass["Neg Pred Value"], 4)
  ),
  F1 = c(
    round(train_cm$byClass["F1"], 4),
    round(test_cm$byClass["F1"], 4)
  ),
  AUC = c(
    round(train_auc, 4),
    round(test_auc, 4)
  )
)

cat("\n===== Performance Metrics Summary =====\n")
print(metrics_table)

write.csv(metrics_table, "results/Model_performance_metrics.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# 9. ROC Curves for Key Genes (Individual Features)
#-------------------------------------------------------------------------------

cat("\n===== Individual Gene ROC Analysis =====\n")

# Key genes identified from ML
key_genes <- c("MAPK8", "SIRT1", "PIK3R1", "KRAS", "MAPK1")

# Check which genes are in the data
available_genes <- key_genes[key_genes %in% colnames(Test_set)]
cat("Available genes for ROC:", paste(available_genes, collapse = ", "), "\n")

if (length(available_genes) > 0) {
  
  gene_roc_data <- list()
  gene_auc_values <- data.frame(Gene = character(), AUC = numeric())
  
  for (gene in available_genes) {
    
    roc_obj <- roc(test_true_num, Test_set[, gene], quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    
    gene_auc_values <- rbind(gene_auc_values, 
                             data.frame(Gene = gene, AUC = round(auc_val, 4)))
    
    # Store ROC curve data
    gene_roc_data[[gene]] <- data.frame(
      Gene = gene,
      TPR = roc_obj$sensitivities,
      FPR = 1 - roc_obj$specificities,
      AUC = round(auc_val, 4)
    )
  }
  
  # Combine ROC data
  all_roc_data <- do.call(rbind, gene_roc_data)
  
  # Color palette
  colors <- c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "grey")
  names(colors) <- available_genes[1:min(5, length(available_genes))]
  
  # Create labels with AUC
  auc_labels <- sapply(available_genes, function(g) {
    auc_val <- gene_auc_values$AUC[gene_auc_values$Gene == g]
    paste0(g, " (AUC = ", auc_val, ")")
  })
  
  p_gene_roc <- ggplot(all_roc_data, aes(x = FPR, y = TPR, color = Gene)) +
    geom_line(linewidth = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = colors, labels = auc_labels) +
    labs(
      title = "ROC Curves - Key Genes",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.7, 0.25),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.text = element_text(size = 10)
    ) +
    coord_fixed() +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1))
  
  ggsave("figures/ROC_key_genes.pdf", p_gene_roc, width = 6, height = 6)
  
  # Save gene AUC table
  write.csv(gene_auc_values, "results/Gene_AUC_values.csv", row.names = FALSE)
  
  cat("\nGene AUC values:\n")
  print(gene_auc_values)
}

#-------------------------------------------------------------------------------
# 10. Save All Results
#-------------------------------------------------------------------------------

evaluation_results <- list(
  best_model = best_model_name,
  train_cm = train_cm,
  test_cm = test_cm,
  train_roc = train_roc,
  test_roc = test_roc,
  train_auc = train_auc,
  test_auc = test_auc,
  metrics = metrics_table
)

saveRDS(evaluation_results, "results/Model_evaluation_results.rds")

cat("\n===== Model Evaluation Completed =====\n")
```
