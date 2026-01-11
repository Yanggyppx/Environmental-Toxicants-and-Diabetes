# Requirements Installation Script.R
# CRAN packages
cran_packages <- c(
  # Data manipulation
  "readxl", "dplyr", "tidyr", "tibble",
  
  # Statistical analysis
  "rms", "car",
  
  # Machine learning
  "caret", "glmnet", "randomForestSRC", "xgboost", "gbm", 
  "e1071", "BART", "MASS", "plsRglm", "mboost",
  
  # Visualization
  "ggplot2", "gridExtra", "patchwork", "RColorBrewer", 
  "colorspace", "ComplexHeatmap", "ggVennDiagram",
  
  # Model interpretation
  "pROC", "shapviz", "fastshap",
  
  # Mixture analysis
  "bkmr", "qgcomp",
  
  # Parallel computing
  "snowfall", "parallel"
)

# Bioconductor packages
bioc_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

cat("All packages installed successfully!\n")